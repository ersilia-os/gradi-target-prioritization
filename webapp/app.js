// app.js — target-prioritization browser (vanilla, no deps).
// Loads per-organism JSON, computes a live weighted composite, filters/sorts,
// renders a windowed OpenTargets-style table, and shows a per-protein detail drawer.
"use strict";

const LS_KEY = "gradi-web2-state";   // distinct from the live site, which shares this origin when served beside it
const PAGE_SIZE = 100; // rows per page

// min-threshold range filters (0–1). "__c" is the live composite.
const RANGES = [];

const state = {
  org: "ec",
  view: "table",          // "table" (one table, expandable column groups) | "compare" | "map"
  expanded: [],           // names of the column groups currently expanded
  cardOpen: [],           // the gene card keeps its own open section, so using the card never changes the table
  mapColorBy: "degradability_consensus",
  mapThr: {},             // map: per colour-by field, the threshold the reader set (else a data-derived default)
  mapLabels: true,        // map: draw cluster labels in categorical colour modes
  weights: {},            // key -> {on, weight}
  orthoTransfer: true,    // include ortholog-transferred knowledge (essentiality track)
  filters: { search: "", ranges: {}, cats: {}, bools: {} },
  shortlist: [],          // starred accessions (basket)
  shortlistOnly: false,   // filter table to the shortlist
  mapSel: null,           // Set of accessions box-selected on the map (or null)
  sort: { key: "__c", dir: -1 },          // the weighted score is on by default (falls back to degradability)
};

const cache = {};         // org -> parsed json
let DATA = null;          // current org payload
let AVAIL = {};           // component key -> bool (has data)
let filtered = [];        // current filtered+sorted rows
let page = 0;   // current page index (0-based)
let weightsDirty = true;  // recompute the composite only when weights actually change
let _lastRow = null;      // row currently shown in the drawer (for live re-render)
let _kbIdx = -1;          // keyboard-focused row index within the current page
function setKbFocus(idx) {
  const rows = [...document.querySelectorAll("#tbody tr")];
  if (!rows.length) return;
  _kbIdx = Math.max(0, Math.min(idx, rows.length - 1));
  rows.forEach((r, i) => r.classList.toggle("kbfocus", i === _kbIdx));
  rows[_kbIdx].scrollIntoView({ block: "nearest" });
}

// Saved or shared state is untrusted: an old link or another version of the site can carry keys that no
// longer exist, and applying them silently filters or colours the table by something that is not on screen.
function cleanWeights(w) {
  const out = {};
  if (!w || typeof w !== "object") return out;
  for (const c of COMPONENTS) {
    const v = w[c.key];
    if (v && typeof v === "object") out[c.key] = { on: !!v.on, weight: isNum(v.weight) ? v.weight : 0 };
  }
  return out;
}
function cleanRanges(r) {
  const out = {};
  if (!r || typeof r !== "object") return out;
  for (const f of RANGE_FILTERS) {
    const v = r[f.key];
    if (Array.isArray(v) && v.length === 2 && v.every(isNum)) out[f.key] = v;
  }
  return out;
}
function cleanCats(c) {
  const out = {};
  if (!c || typeof c !== "object") return out;
  for (const [k, v] of Object.entries(c)) if (Array.isArray(v) && v.length) out[k] = new Set(v.filter((x) => typeof x === "string" || isNum(x)));
  return out;
}
function cleanBools(b) {
  const out = {};
  if (!b || typeof b !== "object") return out;
  for (const [k, v] of Object.entries(b)) if (v === "any" || v === "yes" || v === "no") out[k] = v;
  return out;
}
const cleanMapColor = (k) => (MAP_COLORS.some((m) => m.key === k) ? k : "degradability_consensus");

// ---------- persistence ----------------------------------------------------
function save() {
  const s = {
    org: state.org, view: state.view, mapColorBy: state.mapColorBy,
    mapThr: state.mapThr, mapLabels: state.mapLabels,
    weights: state.weights, orthoTransfer: state.orthoTransfer,
    filters: { search: state.filters.search, ranges: state.filters.ranges,
      cats: Object.fromEntries(Object.entries(state.filters.cats).map(([k, v]) => [k, [...v]])),
      bools: state.filters.bools },
    shortlist: state.shortlist, shortlistOnly: state.shortlistOnly,
    expanded: state.expanded,
    sort: state.sort,
  };
  try { localStorage.setItem(LS_KEY, JSON.stringify(s)); } catch (e) {}
  syncHash();
}
function load() {
  try {
    const s = JSON.parse(localStorage.getItem(LS_KEY) || "null");
    if (!s) return;
    state.org = ORGANISM_META[s.org] ? s.org : "ec";
    state.view = (s.view === "map" || s.view === "compare") ? s.view : "table";
    state.expanded = Array.isArray(s.expanded) ? s.expanded.filter((g) => GROUP_ORDER.includes(g)).slice(0, 1) : [];
    state.mapColorBy = cleanMapColor(s.mapColorBy);
    state.mapThr = (s.mapThr && typeof s.mapThr === "object") ? s.mapThr : {};
    if (typeof s.mapLabels === "boolean") state.mapLabels = s.mapLabels;
    state.weights = cleanWeights(s.weights);
    if (typeof s.orthoTransfer === "boolean") state.orthoTransfer = s.orthoTransfer;
    if (s.filters) {
      state.filters.search = s.filters.search || "";
      state.filters.ranges = cleanRanges(s.filters.ranges);
      state.filters.cats = cleanCats(s.filters.cats);
      state.filters.bools = cleanBools(s.filters.bools);
    }
    state.shortlist = Array.isArray(s.shortlist) ? s.shortlist : [];
    state.shortlistOnly = !!s.shortlistOnly;
    state.sort = s.sort || state.sort;
  } catch (e) {}
}

// ---------- shareable deep-link (URL hash) ---------------------------------
let _applyingHash = false;
function serializeState() {
  const arr = (o) => Object.fromEntries(Object.entries(o).filter(([, v]) => v && v.length));
  return {
    org: state.org, view: state.view, ex: state.expanded, map: state.mapColorBy,
    w: state.weights,
    rg: Object.keys(state.filters.ranges).length ? state.filters.ranges : undefined,
    q: state.filters.search || undefined,
    cats: arr(Object.fromEntries(Object.entries(state.filters.cats).map(([k, v]) => [k, [...v]]))),
    bools: Object.fromEntries(Object.entries(state.filters.bools).filter(([, m]) => m && m !== "any")),
    sl: state.shortlist.length ? state.shortlist : undefined,
    so: state.shortlistOnly || undefined,
    sel: state._sel || undefined,
  };
}
function syncHash() {
  if (_applyingHash) return;
  try { history.replaceState(null, "", "#p=" + btoa(encodeURIComponent(JSON.stringify(serializeState())))); } catch (e) {}
}
function applyHash() {
  const m = location.hash.match(/[#&]p=([^&]+)/);
  if (!m) return false;
  try {
    const p = JSON.parse(decodeURIComponent(atob(m[1])));
    _applyingHash = true;
    if (p.org && ORGANISM_META[p.org]) state.org = p.org;
    if (p.view === "map" || p.view === "compare") state.view = p.view; else state.view = "table";
    if (Array.isArray(p.ex)) state.expanded = p.ex.filter((g) => GROUP_ORDER.includes(g)).slice(0, 1);
    if (typeof p.xfer === "boolean") state.orthoTransfer = p.xfer;
    if (p.map) state.mapColorBy = cleanMapColor(p.map);
    if (p.w) state.weights = cleanWeights(p.w);
    state.filters.ranges = cleanRanges(p.rg);
    state.filters.search = p.q || "";
    state.filters.cats = cleanCats(p.cats);
    state.filters.bools = cleanBools(p.bools);
    state.shortlist = p.sl || [];
    state.shortlistOnly = !!p.so;
    if (p.sel) state._sel = p.sel;
    _applyingHash = false;
    return true;
  } catch (e) { _applyingHash = false; return false; }
}

// ---------- helpers --------------------------------------------------------
const $ = (id) => document.getElementById(id);
const isNum = (v) => typeof v === "number" && !Number.isNaN(v);

function fmt(type, v) {
  if (v === null || v === undefined || v === "") return "–";
  switch (type) {
    case "score": return isNum(v) ? v.toFixed(2) : String(v);
    case "pct": return isNum(v) ? Math.round(v * 100) + "%" : String(v);
    case "plddt": return isNum(v) ? v.toFixed(0) : String(v);
    case "num": return isNum(v) ? (+v).toFixed(2) : String(v);
    case "int": return isNum(v) ? String(Math.round(v)) : String(v);
    default: return String(v);
  }
}
function heatStyle(v, color) {
  if (!isNum(v)) return "";
  const c = color || "var(--brand)";
  const amt = Math.max(0, Math.min(100, Math.round(12 + v * 78)));
  return `background:color-mix(in srgb,${c} ${amt}%,#fff);`
       + (amt >= 58 ? "color:#fff;" : "color:var(--ink);")
       + (v === 0 ? "box-shadow:inset 0 0 0 1px var(--border);color:var(--faint);" : "");  // true zero ≠ no data
}
function badgeHTML(type, v, key) {
  if (v === null || v === undefined || v === "") return "–";
  const pfx = type === "sel" ? "sel-" : type === "pop" ? "pop-" : type === "loc" ? "loc-" : "tier-";
  const full = String(v).replace(/_/g, " ");
  const abbr = (key && CAT_ABBREV[key] && CAT_ABBREV[key][v]) || full.slice(0, 2);
  return `<span class="badge ${pfx}${v}" title="${full}">${abbr}</span>`;
}
function boolHTML(v) {
  if (v === true) return '<span class="yes">✓</span>';
  if (v === false) return '<span class="no">·</span>';
  return '<span class="no">–</span>';
}

// ---------- composite ------------------------------------------------------
// the weighted-score column only exists while at least one slider is on
function showWeighted() { return enabledComponents().length > 0; }
function enabledComponents() {
  return COMPONENTS.filter((c) => AVAIL[c.key] && state.weights[c.key] && state.weights[c.key].on
    && state.weights[c.key].weight > 0);
}
function computeAll() {
  const comps = enabledComponents();
  for (const row of DATA.rows) {
    let sum = 0, wsum = 0;
    for (const c of comps) {
      let v = row[c.key];
      if (isNum(v) && c.invert) v = 1 - v;
      if (isNum(v)) { const w = state.weights[c.key].weight; sum += w * v; wsum += w; }
    }
    row.__c = wsum > 0 ? sum / wsum : null;
  }
}
// ---------- distribution / percentile (per org, cached) --------------------
let _pctCache = {};
function percentileOf(key, v) {
  if (!isNum(v)) return null;
  let arr = _pctCache[key];
  if (!arr) { arr = DATA.rows.map((r) => r[key]).filter(isNum).sort((a, b) => a - b); _pctCache[key] = arr; }
  if (!arr.length) return null;
  let lo = 0, hi = arr.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] <= v) lo = m + 1; else hi = m; }
  return Math.round(100 * lo / arr.length);
}
function tieShare(key, v) {
  const arr = _pctCache[key];
  if (!arr || !arr.length) return 0;
  let lo = 0, hi = arr.length;
  while (lo < hi) { const m = (lo + hi) >> 1; if (arr[m] < v) lo = m + 1; else hi = m; }
  let up = lo; while (up < arr.length && arr[up] === v) up++;
  return (up - lo) / arr.length;
}
function ordinal(n) {
  const s = ["th", "st", "nd", "rd"], v = n % 100;
  return n + (s[(v - 20) % 10] || s[v] || s[0]);
}

// (The v1 orthology-transfer toggle is gone: v2 tables carry the transfer inside each axis.)
const TRANSFER_KEYS = new Set();
function stashEssentiality() {}
function applyTransferMode() {}

// ---------- filtering / sorting -------------------------------------------
function passFilters(row, skipKey) {
  const f = state.filters;
  if (state.shortlistOnly && !state.shortlist.includes(row.uniprot_accession)) return false;
  if (state.mapSel && !state.mapSel.has(row.uniprot_accession)) return false;
  if (f.search) {
    const q = f.search.toLowerCase();
    const hay = [row.gene, row.uniprot_accession, row.name].filter(Boolean).join(" ").toLowerCase();
    if (!hay.includes(q)) return false;
  }
  for (const key in f.ranges) {
    if (key === skipKey) continue;
    const [lo, hi] = f.ranges[key], v = row[key];
    if (!isNum(v) || v < lo - 1e-9 || v > hi + 1e-9) return false;
  }
  for (const [key, set] of Object.entries(f.cats)) {
    if (key === skipKey || !set || set.size === 0) continue;
    if (!set.has(row[key])) return false;
  }
  for (const [key, mode] of Object.entries(f.bools)) {
    if (!mode || mode === "any") continue;
    const want = mode === "yes";
    if (Boolean(row[key]) !== want) return false;
  }
  return true;
}
let _sorted = null, _sortSig = "";
function sortedRows() {
  const { dir } = state.sort;
  const key = state.sort.key === "__c" && !showWeighted() ? "degradability_consensus" : state.sort.key;
  const sig = `${state.org}|${key}|${dir}|${DATA.__ver || 0}`;
  if (_sorted && sig === _sortSig) return _sorted;
  _sorted = DATA.rows.slice().sort((a, b) => {
    const x = a[key], y = b[key];
    const xn = x === null || x === undefined || x === "", yn = y === null || y === undefined || y === "";
    if (xn && yn) return 0;
    if (xn) return 1;      // nulls always last
    if (yn) return -1;
    if (isNum(x) && isNum(y)) return (x - y) * dir;
    return String(x).localeCompare(String(y)) * dir;
  });
  _sortSig = sig;
  return _sorted;
}
let _enter = false;
function recompute(opts) {
  const quiet = !!(opts && opts.quiet);
  _enter = !quiet;
  if (weightsDirty) { computeAll(); weightsDirty = false; DATA.__ver = (DATA.__ver || 0) + 1; }
  filtered = sortedRows().filter((r) => passFilters(r));
  page = 0;
  if (!quiet) $("tablewrap").scrollTop = 0;
  $("nShown").textContent = filtered.length.toLocaleString();
  $("sbDoneN").textContent = filtered.length.toLocaleString(); $("sbShown").textContent = filtered.length.toLocaleString(); $("sbTotal").textContent = DATA.rows.length.toLocaleString();
  updateAccCounts();
  $("nTotal").textContent = DATA.rows.length.toLocaleString();
  renderActiveView();
  refreshAll();
}
function renderActiveView() {
  const isMap = state.view === "map", isCmp = state.view === "compare";
  $("tablewrap").hidden = isMap || isCmp;
  $("mapwrap").hidden = !isMap;
  $("cmpwrap").hidden = !isCmp;
  $("mapctrl").hidden = !isMap;
  $("pager").classList.toggle("ghost", isMap || isCmp);
  if (isMap) { drawMap(); const pn = $("provNote"); if (pn) pn.hidden = true; }
  else if (isCmp) { renderCompare(); const pn = $("provNote"); if (pn) pn.hidden = true; }
  else { renderThead(); renderPage(); }
  updateMapSelChip();
  setTopbarOffset();
}
// side-by-side comparison of the shortlisted targets (metrics as rows, targets as columns)
const CMP_KEYS = ["__c", "length", "degradability_consensus", "adep4_prob", "onc212_prob",
  "essentiality_consensus", "essentiality_evidence", "deg_essential_any", "ligands_consensus", "ligands_evidence", "n_drugs",
  "pockets_consensus", "localization", "cytoplasmic_fraction", "has_human_ortholog", "bacterial_panel_orthologs",
  "studiedness_consensus", "cog"];
function renderCompare() {
  const wrap = $("cmpwrap");
  const rows = state.shortlist.map((a) => DATA.rows.find((r) => r.uniprot_accession === a)).filter(Boolean);
  if (!rows.length) {
    wrap.innerHTML = `<div class="cmpempty">Star proteins to compare them side by side.<span>Use the star beside a protein in the table or in its card.</span></div>`;
    return;
  }
  const byKey = Object.fromEntries(TABLE_COLUMNS.map((c) => [c.key, c]));
  const label = (k) => k === "__c" ? "Weighted (yours)" : !byKey[k] ? k
    : (["Rank", "Evidence", "Class", "Measured"].includes(byKey[k].label) ? `${byKey[k].group} ${byKey[k].label.toLowerCase()}` : byKey[k].label);
  const keys = CMP_KEYS.filter((k) => (k === "__c" && showWeighted()) || (k !== "__c" && byKey[k] && DATA.columns.includes(k)));
  const cell = (r, k) => {
    if (k === "__c") return `<span class="heat" style="${heatStyle(isNum(r.__c) ? r.__c : null, AXIS_COLORS.composite)}">${isNum(r.__c) ? r.__c.toFixed(2) : "–"}</span>`;
    return cellHTML(r, byKey[k]);
  };
  let h = `<table class="cmptable"><thead><tr><th class="cmlh"></th>`
    + rows.map((r) => `<th><div class="gene">${esc(r.gene || r.uniprot_accession)}</div><div class="acc">${esc(r.uniprot_accession)}</div></th>`).join("")
    + `</tr></thead><tbody>`;
  for (const k of keys) {
    h += `<tr><td class="cml" style="--ac:${colColor(k)}">${label(k)}</td>`
      + rows.map((r) => `<td class="rc">${cell(r, k)}</td>`).join("") + `</tr>`;
  }
  wrap.innerHTML = h + `</tbody></table>`;
}
// segmented controls slide a thumb behind the selected button
function syncSeg(host) {
  if (!host) return;
  const bs = [...host.querySelectorAll("button")];
  const i = Math.max(0, bs.findIndex((b) => b.getAttribute("aria-selected") === "true"));
  host.style.setProperty("--i", i); host.style.setProperty("--n", bs.length);
}
function buildViewSwitch() {
  const host = $("viewSwitch");
  const tabs = [{ view: "table", label: "Table" }, { view: "compare", label: "⇄ Compare" }, { view: "map", label: "◵ Map" }];
  host.innerHTML = tabs.map((t) =>
    `<button role="tab" data-view="${t.view}" aria-selected="${state.view === t.view}">${t.label}</button>`).join("");
  host.querySelectorAll("button").forEach((b) => b.onclick = () => applyView(b.dataset.view));
  syncSeg(host);
}
function applyView(key) {
  state.view = key;
  $("viewSwitch").querySelectorAll("button").forEach((x) =>
    x.setAttribute("aria-selected", x.dataset.view === key));
  syncSeg($("viewSwitch"));
  buildFilterBar();
  save();
  recompute();
}
// ---------- column groups: one table, each group expands / collapses ------------
// One group open at a time: opening another closes the current one; clicking the open group closes it.
function toggleGroup(name) {
  state.expanded = state.expanded[0] === name ? [] : [name];
  save(); renderActiveView();
}
// A column with no value at all for this organism is dropped (e.g. measured essentiality on K. pneumoniae).
function colHasData(k) {
  const m = DATA.__has || (DATA.__has = {});
  if (m[k] === undefined) m[k] = DATA.rows.some((r) => r[k] !== null && r[k] !== undefined);
  return m[k];
}
function visibleColumns() {
  const out = [];
  for (const g of GROUP_ORDER) {
    const open = state.expanded.includes(g);
    for (const c of TABLE_COLUMNS) {
      if (c.group !== g || !DATA.columns.includes(c.key) || !colHasData(c.key)) continue;
      if (open || c.summary) out.push(c);
    }
  }
  return out.concat(TABLE_COLUMNS.filter((c) => c.fixed && DATA.columns.includes(c.key)));   // Length last
}
function tipAttr(text) {
  return text ? ` data-tip="${String(text).replace(/"/g, "&quot;")}"` : "";
}
const ROTATE_TYPES = new Set(["score", "pct", "plddt", "num", "int", "bool", "binary01"]);
function colWidths() {
  const w = window.innerWidth;
  if (w <= 700) return { rank: 34, gene: 150, data: 100, wide: 190 };
  if (w <= 900) return { rank: 40, gene: 190, data: 100, wide: 210 };
  if (w <= 1440) return { rank: 44, gene: 210, data: 98, wide: 210, k: 0.94 };   // laptop: about a quarter more columns in view
  return { rank: 44, gene: 240, data: 108, wide: 230 };
}
// a collapsed group is one column wide and its band carries the group name, so that column never goes below 116
function colW(c, W) {
  const base = c.w || (c.type === "text" ? W.wide : W.data);
  if (c.summary) return Math.max(base, 116);
  return Math.round(base * (c.w && W.k && c.key !== "length" ? W.k : 1));
}
function renderThead() {
  const cols = visibleColumns();
  const arrow = (k) => `<span class="arrow">${state.sort.key === k ? (state.sort.dir < 0 ? "↓" : "↑") : ""}</span>`;
  const th = (key, label, { nosort, tip, cls: extra } = {}) => {
    const sortAttr = nosort ? "" : ` data-sort="${key}"`;
    const ariaSort = (!nosort && key === state.sort.key) ? ` aria-sort="${state.sort.dir < 0 ? "descending" : "ascending"}"` : "";
    return `<th${sortAttr}${ariaSort} class="${[nosort ? "nosort" : "", extra || ""].join(" ").trim()}"${tipAttr(tip)}><span class="hc"><span class="hl">${label}</span>${nosort ? "" : arrow(key)}</span></th>`;
  };
  let h = th("", "#", { nosort: true, tip: LEADING_COL_DESC.rank, cls: "colrank" })
    + th("gene", "Protein", { tip: LEADING_COL_DESC.name, cls: "colgene" })
    + (showWeighted() ? th("__c", "Weighted", { tip: LEADING_COL_DESC.__c, cls: "ar" }) : "");
  // a group's headline score needs no label while the group is collapsed: the band above already names it
  const headLabel = (c) => !c.key.endsWith("_consensus") ? c.label : state.expanded.includes(c.group) ? "Score" : "";
  for (const c of cols) h += th(c.key, headLabel(c), { tip: c.desc, cls: alignOf(c) });

  // group band: merge consecutive columns of one group; the untitled run (length, weighted) stays blank
  const runs = [];
  for (const col of [...(showWeighted() ? [{ key: "__c" }] : []), ...cols]) {
    const g = columnGroup(col);
    const last = runs[runs.length - 1];
    if (last && last.group === g) last.span++;
    else runs.push({ group: g, span: 1 });
  }
  let band = `<th class="gband gband-empty gband-corner" colspan="2"></th>`;
  for (const r of runs) {
    if (!r.group) { band += `<th class="gband gband-empty" colspan="${r.span}"></th>`; continue; }
    const open = state.expanded.includes(r.group);
    // a collapsed group is one column wide: long names take their short form there
    const label = (!open && r.group.length > 14 && GROUP_DISPLAY[r.group]) || r.group;
    band += `<th class="gband gtoggle${open ? " open" : ""}" data-group="${r.group}" colspan="${r.span}" style="--gc:${groupColor(r.group)}"`
      + tipAttr(GROUP_NOTES[r.group]) + ` role="button" aria-expanded="${open}"><span class="gl"><i class="chev"></i>${label}</span></th>`;
  }

  // explicit column widths (table-layout: fixed reads only the first row, which has a colspan)
  const W = colWidths();
  const widths = [W.rank, W.gene, ...(showWeighted() ? [W.data] : []), ...cols.map((c) => (c.key === "length" && window.innerWidth <= 700) ? W.data : colW(c, W))];
  document.documentElement.style.setProperty("--rankw", W.rank + "px");
  const tbl = $("table");
  let cg = tbl.querySelector("colgroup");
  if (!cg) { cg = document.createElement("colgroup"); tbl.insertBefore(cg, tbl.firstChild); }
  cg.innerHTML = widths.map((w) => `<col style="width:${w}px">`).join("");
  tbl.style.width = widths.reduce((s, w) => s + w, 0) + "px";
  const thead = document.querySelector("#table thead");
  thead.innerHTML = `<tr class="bandrow">${band}</tr><tr>${h}</tr>`;
  thead.querySelectorAll("th.gtoggle").forEach((el) => { el.onclick = () => toggleGroup(el.dataset.group); });
  thead.querySelectorAll("th[data-sort]").forEach((el) => {
    el.onclick = () => {
      const k = el.dataset.sort;
      if (state.sort.key === k) state.sort.dir *= -1;
      else state.sort = { key: k, dir: (k === "gene") ? 1 : -1 };
      save(); recompute();
    };
  });
}
// numbers align right, short categorical marks centre, words left
function alignOf(c) {
  if (["score", "pct", "plddt", "num", "int"].includes(c.type)) return "ar";
  if (["tier", "bool", "binary01", "class"].includes(c.type)) return "ac";
  return "al";
}
function tierType(key) {
  return key === "localization" ? "loc" : "tier";
}
const LOC_LABEL = { cytoplasm: "Cytoplasm", cytoplasmic_membrane: "Inner membrane", periplasm: "Periplasm", outer_membrane: "Outer membrane", extracellular: "Extracellular" };
function classBadgeHTML(v) {
  if (!v || v === "none") return '<span class="na">–</span>';
  const m = FC_BY_ID[v];
  return `<span class="cog" title="${esc(m ? m.label : String(v))}">${esc(v)}</span>`;
}
function cellHTML(row, c) {
  const v = row[c.key];
  const NA = '<span class="na">–</span>';
  if (v === null || v === undefined || v === "") return NA;
  if (c.type === "bool" || c.type === "binary01") {
    const yes = v === true || v === 1 || (typeof v === "number" && v >= 0.5);
    if (c.key === "has_human_ortholog") return yes ? '<span class="txt">yes</span>' : '<span class="na">no</span>';
    return yes ? '<span class="yes">✓</span>' : '<span class="na">–</span>';
  }
  if (c.type === "class") return classBadgeHTML(v);
  if (c.type === "loc") return `<span class="txt" title="${esc(LOC_LABEL[v] || v)}">${esc(LOC_LABEL[v] || String(v).replace(/_/g, " "))}</span>`;
  if (c.type === "tier") return `<span class="ev e${v}" title="Evidence level ${v} of 3">${v}</span>`;
  if (c.type === "score" && c.heat && isNum(v))
    return `<span class="sc${v === 0 ? " zero" : ""}" style="--c:${colColor(c.key)};--v:${Math.round(Math.max(0, Math.min(1, v)) * 100)}%"><b>${v.toFixed(2)}</b></span>`;
  if (c.type === "text") {
    const full = String(v).replace(/;/g, "; ");
    return `<span class="celltext" title="${full.replace(/"/g, "&quot;")}">${esc(full)}</span>`;
  }
  return `<span class="numcell">${fmt(c.type === "int" && c.key === "length" ? "int" : c.type, v)}</span>`;
}
function pageCount() { return Math.max(1, Math.ceil(filtered.length / PAGE_SIZE)); }
let _renderTok = 0;
const FIRST_ROWS = 28;   // enough to fill the screen; the rest follows next frame
function renderPage() {
  const cols = visibleColumns();
  if (page >= pageCount()) page = pageCount() - 1;
  if (page < 0) page = 0;
  const start = page * PAGE_SIZE;
  const end = Math.min(start + PAGE_SIZE, filtered.length);
  const frag = document.createDocumentFragment();
  if (!filtered.length) {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td class="emptyrow" colspan="${cols.length + 2 + (showWeighted() ? 1 : 0)}">No proteins match. <button class="linkbtn" id="emptyClear">Clear filters</button></td>`;
    frag.appendChild(tr);
  }
  let html = "";
  for (let i = start; i < end; i++) {
    const row = filtered[i];
    const cc = isNum(row.__c) ? row.__c : null;
    const pdesc = row.name || "";
    html += `<tr data-acc="${row.uniprot_accession}"${row.uniprot_accession === state._sel ? ' class="sel"' : ""}>`
      + `<td class="colrank rank">${i + 1}</td>`
      + `<td class="colgene"${pdesc ? ` title="${pdesc.replace(/"/g, "&quot;")}"` : ""}>${starHTML(row.uniprot_accession)}`
      + `<div class="genehead">${row.gene ? `<span class="gene">${esc(row.gene)}</span>` : ""}`
      + `<span class="acc${row.gene ? "" : " accsolo"}">${row.uniprot_accession}</span></div>`
      + (pdesc ? `<div class="pdesc">${esc(pdesc)}</div>` : "") + `</td>`
      + (showWeighted() ? `<td class="rc ar">${cellHTML({ __c: cc }, { key: "__c", type: "score", heat: true })}</td>` : "");
    for (const col of cols) html += `<td class="rc ${alignOf(col)}">${cellHTML(row, col)}</td>`;
    html += "</tr>";
  }
  const tb = $("tbody");
  const tok = ++_renderTok;
  if (!filtered.length) tb.innerHTML = frag.firstChild ? frag.firstChild.outerHTML : "";
  else {
    let cut = 0;
    for (let k = 0; k < FIRST_ROWS; k++) cut = html.indexOf("</tr>", cut) + 5;   // end of the FIRST_ROWS-th row
    if (end - start <= FIRST_ROWS) tb.innerHTML = html;
    else {
      tb.innerHTML = html.slice(0, cut);
      requestAnimationFrame(() => { if (tok === _renderTok) tb.insertAdjacentHTML("beforeend", html.slice(cut)); });
    }
  }
  if (!tb._wired) {
    tb._wired = true;
    tb.addEventListener("click", (e) => {
      const tr = e.target.closest("tr[data-acc]");
      if (!tr) return;
      const acc = tr.dataset.acc;
      if (e.target.closest(".star")) { e.stopPropagation(); toggleStar(acc); return; }
      const row = (DATA.__byAcc || (DATA.__byAcc = new Map(DATA.rows.map((r) => [r.uniprot_accession, r])))).get(acc);
      if (row) openDrawer(row);
    });
  }
  const ec = $("emptyClear"); if (ec) ec.onclick = () => applyPreset("reset");
  _kbIdx = -1;
  $("table").hidden = false;
  $("loading").hidden = true;
  tb.classList.remove("enter");
  if (_enter) { void tb.offsetWidth; tb.classList.add("enter"); _enter = false; }
  renderPager();
  setTopbarOffset();
}
// ---------- methods / provenance modal -------------------------------------
function openMethods() {
  const m = $("methodsModal");
  const refsHTML = (refs) => !refs || !refs.length ? "" :
    `<ul class="mrefs">` + refs.map(([txt, url]) => `<li>${url
      ? `<a href="${url}" target="_blank" rel="noopener">${txt} ↗</a>` : txt}</li>`).join("") + `</ul>`;
  m.innerHTML = `<div class="modalhead"><h2>How to read each axis</h2><button class="close" id="methodsClose" aria-label="Close">×</button></div>`
    + `<div class="modalbody">`
    + METHODS.map((s) => `<section class="mrow" style="--pc:${s.color}"><h3>${s.title}</h3><p>${s.body}</p>${refsHTML(s.refs)}</section>`).join("")
    + `<p class="mfoot">Full axis specifications live in the project <code>docs/&lt;axis&gt;.md</code>. This is a preview built from the v2 tables.</p>`
    + `</div>`;
  m.hidden = false; $("methodsScrim").hidden = false;
  $("methodsClose").onclick = closeMethods;
}
function closeMethods() { $("methodsModal").hidden = true; $("methodsScrim").hidden = true; }
// ---------- shortlist / basket ---------------------------------------------
function isStarred(acc) { return state.shortlist.includes(acc); }
function starHTML(acc) {
  return `<button class="star${isStarred(acc) ? " on" : ""}" data-acc="${acc}" title="Star to compare" aria-label="Star">★</button>`;
}
function toggleStar(acc) {
  const i = state.shortlist.indexOf(acc);
  if (i >= 0) state.shortlist.splice(i, 1); else state.shortlist.push(acc);
  save(); updateShortlistBtn();
  if (state.shortlistOnly || state.view === "compare") recompute();
  else document.querySelectorAll(`.star[data-acc="${acc}"]`).forEach((el) => el.classList.toggle("on", isStarred(acc)));
}
function updateShortlistBtn() {
  const b = $("shortlistBtn");
  if (!b) return;
  b.innerHTML = `★<span class="lbl"> Starred</span> (${state.shortlist.length})`;
  b.setAttribute("aria-pressed", state.shortlistOnly ? "true" : "false");
  b.classList.toggle("on", state.shortlistOnly);
}
function updateMapSelChip() {
  const c = $("mapSelChip");
  if (!c) return;
  if (state.mapSel && state.mapSel.size) {
    c.hidden = false;
    c.innerHTML = `${state.mapSel.size.toLocaleString()} on map <button id="mapSelClear" title="Clear map selection">✕</button>`;
    $("mapSelClear").onclick = (e) => { e.stopPropagation(); state.mapSel = null; updateMapSelChip(); recompute(); };
  } else { c.hidden = true; }
}
function gotoPage(p) {
  const n = pageCount();
  page = Math.max(0, Math.min(p, n - 1));
  renderPage();
  $("tablewrap").scrollTop = 0;
}
function renderPager() {
  const foot = $("pagerFoot");
  if (foot) {
    foot.hidden = filtered.length <= PAGE_SIZE || state.view !== "table";
    foot.innerHTML = `<button class="pbtn" data-p="prev" ${page === 0 ? "disabled" : ""}>‹ Previous</button><span class="pinfo">page <b>${page + 1}</b> of ${pageCount()}</span><button class="pbtn" data-p="next" ${page >= pageCount() - 1 ? "disabled" : ""}>Next ›</button>`;
    foot.querySelectorAll("button").forEach((b) => b.onclick = () => gotoPage(b.dataset.p === "prev" ? page - 1 : page + 1));
  }
  const host = $("pager");
  const n = pageCount();
  host.classList.toggle("ghost", state.view !== "table");
  if (!filtered.length) { host.innerHTML = `<span class="pinfo">no matches</span>`; return; }
  host.innerHTML =
      `<button class="pbtn" data-p="prev" title="Previous page" ${page === 0 ? "disabled" : ""}>‹</button>`
    + `<span class="pinfo">page <b>${(page + 1).toLocaleString()}</b> / ${n.toLocaleString()}</span>`
    + `<button class="pbtn" data-p="next" title="Next page" ${page >= n - 1 ? "disabled" : ""}>›</button>`;
  host.querySelectorAll("button[data-p]").forEach((b) => {
    b.onclick = () => gotoPage(b.dataset.p === "prev" ? page - 1 : page + 1);
  });
}
// measure the sticky top bar so column headers can stick right beneath it
function setTopbarOffset() {
  const tb = document.getElementById("topbar");
  if (tb) document.documentElement.style.setProperty("--tbh", tb.offsetHeight + "px");
}

// ---------- KPI tiles ------------------------------------------------------

// ---------- detail drawer --------------------------------------------------
// ---- gene-card render helpers --------------------------------------------
function esc(s) {
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
}
// composite rank among all rows under the current weighting
function compRank(row) {
  const c = row.__c;
  if (!isNum(c)) return null;
  let r = 1;
  for (const x of DATA.rows) if (isNum(x.__c) && x.__c > c) r++;
  return r;
}
// ---- gene card: the same column groups as the table, same open/closed state --------------------
function groupPanelsHTML(row) {
  let h = "";
  for (const g of GROUP_ORDER) {
    const cols = TABLE_COLUMNS.filter((c) => c.group === g && DATA.columns.includes(c.key) && colHasData(c.key));
    if (!cols.length) continue;
    const open = state.cardOpen.includes(g);
    const head = cols.find((c) => c.summary);
    const evid = cols.find((c) => c.type === "tier" && /_evidence$/.test(c.key));
    // three fixed slots (value, percentile, evidence) so the right edge of every group header lines up
    let mv = "", mp = "", me = "";
    if (head && row[head.key] !== null && row[head.key] !== undefined) {
      mv = cellHTML(row, { ...head, heat: false });
      if (head.type === "score" && isNum(row[head.key])) {
        let pc = percentileOf(head.key, row[head.key]);
        if (pc !== null && tieShare(head.key, row[head.key]) > 0.25) pc = null;
        if (pc !== null) mp = `<span class="pctnote">${ordinal(pc)} percentile</span>`;
      }
    }
    if (evid && row[evid.key] != null) me = `<span class="evlab">evidence</span>${cellHTML(row, evid)}`;
    const meta = `<span class="mv">${mv}</span><span class="mp">${mp}</span><span class="me">${me}</span>`;
    let body = "";
    if (open) {
      if (g === "Essentiality" && !(isNum(row.deg_n_datasets) && row.deg_n_datasets > 0))
        body += `<p class="gcnote">Not measured on this strain; every value here is a prediction.</p>`;
      body += `<div class="gcrows">` + cols.filter((c) => c !== head && c !== evid).map((c) => {   // the header line already shows these two
        const empty = row[c.key] === null || row[c.key] === undefined || row[c.key] === "";
        return `<div class="gcrow"${tipAttr(c.desc)}><span class="gck">${esc(c.label)}</span><span class="gcv">${empty ? '<span class="na">–</span>' : cellHTML(row, c)}</span></div>`;
      }).join("") + `</div>`;
      if (g === "Function" && row.cog_categories) body += `<div class="gcrow"><span class="gck">COG categories</span><span class="gcv">`
        + String(row.cog_categories).split(";").filter(Boolean).map((l) => `<span class="cog" title="${esc((FC_BY_ID[l] || {}).label || "")}">${esc(l)}</span>`).join(" ") + `</span></div>`;
    }
    h += `<section class="gcard${open ? " open" : ""}" style="--pc:${groupColor(g)}"><button class="gchead" data-group="${g}" aria-expanded="${open}">`
      + `<i class="chev"></i><span class="gctitle">${g}</span><span class="gcmeta">${meta}</span></button>${body}</section>`;
  }
  return h;
}
function has_(k) { return DATA.columns.includes(k); }
function renderGroupPanels(row) {
  const host = $("groupPanels");
  if (!host) return;
  host.innerHTML = groupPanelsHTML(row);
  host.querySelectorAll(".gchead").forEach((b) => b.onclick = () => {
    state.cardOpen = state.cardOpen[0] === b.dataset.group ? [] : [b.dataset.group];
    renderGroupPanels(row);   // only the card changes
  });
}
// ---- 3D structure viewer (3Dmol.js, lazy) --------------------------------
let _3dmolP = null;       // promise: 3Dmol library loaded
let _mol = null;          // { el, viewer } — one persistent WebGL viewer, reused
const pocketCache = {};   // org -> {acc: [resids]}
function load3Dmol() {
  if (window.$3Dmol) return Promise.resolve();
  if (_3dmolP) return _3dmolP;
  _3dmolP = new Promise((res, rej) => {
    const s = document.createElement("script");
    s.src = "vendor/3Dmol-min.js"; s.async = true;
    s.onload = () => res(); s.onerror = () => { _3dmolP = null; rej(new Error("3Dmol load failed")); };
    document.head.appendChild(s);
  });
  return _3dmolP;
}
// v2 ships pocket SCORES, not pocket residues, so there is no overlay to load.
async function loadPockets(org) { return (pocketCache[org] = {}); }
function ensureMolViewer(slot) {
  if (!window.$3Dmol) return null;
  if (!_mol) {
    const el = document.createElement("div");
    el.className = "mol3dhost";
    slot.appendChild(el);
    _mol = { el, viewer: $3Dmol.createViewer(el, { backgroundColor: "0xf4f4f8", antialias: true }) };
  } else if (_mol.el.parentNode !== slot) {
    slot.appendChild(_mol.el);
  }
  return _mol;
}
function plddtAtomColor(atom) {
  const b = atom.b;
  return b >= 90 ? 0x0053D6 : b >= 70 ? 0x65CBF3 : b >= 50 ? 0xFFDB13 : 0xFF7D45;
}
async function fetchAFModel(acc) {
  for (const v of ["v6", "v4"]) {
    try {
      const r = await fetch(`https://alphafold.ebi.ac.uk/files/AF-${acc}-F1-model_${v}.pdb`);
      if (r.ok) return await r.text();
    } catch (e) {}
  }
  return null;
}
async function renderStructure(row) {
  const slot = document.getElementById("pdbslot");
  if (!slot) return;
  const acc = row.uniprot_accession;
  const msg = (html) => { if (document.getElementById("pdbslot") === slot) slot.innerHTML = `<div class="pdbmsg">${html}</div>`; };
  if (!isNum(row.af_plddt)) {
    msg(`No AlphaFold DB model for this protein (it may be an ESMFold model or too short or long).`);
    const cap = document.getElementById("pdbcap"); if (cap) cap.textContent = "";
    return;
  }
  try { await load3Dmol(); } catch (e) { msg("3D viewer could not be loaded."); return; }
  const [pdb, pockets] = await Promise.all([fetchAFModel(acc), loadPockets(state.org)]);
  if (state._sel !== acc) return;                    // user moved on while awaiting
  if (!pdb || !window.$3Dmol) {
    msg(`Structure unavailable · <a href="https://alphafold.ebi.ac.uk/entry/${acc}" target="_blank" rel="noopener">open in AlphaFold ↗</a>`);
    return;
  }
  slot.innerHTML = "";
  // wait for the slot to have a real layout size before creating/resizing the
  // WebGL viewer, otherwise 3Dmol renders into a zero-size framebuffer (warnings)
  const paint = () => {
    if (document.getElementById("pdbslot") !== slot || state._sel !== acc) return;
    if (!slot.clientWidth || !slot.clientHeight) { requestAnimationFrame(paint); return; }
    const m = ensureMolViewer(slot);
    if (!m) { msg("3D viewer unavailable."); return; }
    const v = m.viewer;
    v.resize();
    v.clear();
    v.addModel(pdb, "pdb");
    v.setStyle({}, { cartoon: { colorfunc: plddtAtomColor } });
    const resids = pockets && pockets[acc];
    if (resids && resids.length) {
      const sel = { resi: resids };
      v.setStyle(sel, { cartoon: { color: 0x2FBF71 } });
      v.addStyle(sel, { stick: { color: 0x1E7A45, radius: 0.18 } });
      try { v.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.62, color: 0x2FBF71 }, sel); } catch (e) {}
    }
    v.zoomTo();      // always frame the WHOLE protein…
    v.zoom(0.9);     // …with a small margin so it never crops
    v.render();
    const cap = document.getElementById("pdbcap");
    if (cap) cap.innerHTML = (resids && resids.length)
      ? `AlphaFold model · cartoon coloured by pLDDT. Drag to rotate, scroll to zoom.`
      : `AlphaFold model · cartoon coloured by pLDDT. Drag to rotate, scroll to zoom.`;
  };
  requestAnimationFrame(paint);
}
function openDrawer(row) {
  state._sel = row.uniprot_accession;
  _lastRow = row;
  document.querySelectorAll("#tbody tr").forEach((tr) =>
    tr.classList.toggle("sel", tr.dataset.acc === row.uniprot_accession));
  const fcMeta = row.cog && row.cog !== "none" ? FC_BY_ID[row.cog] : null;
  const fcPill = fcMeta ? `<span class="fcpill" style="--fc:${fcMeta.color}" title="First COG category">${fcMeta.label}</span>` : "";
  const links = externalLinks(row).map((l) =>
    `<a href="${l.href}" target="_blank" rel="noopener">${l.icon ? `<span class="li">${l.icon}</span>` : ""}${l.label}</a>`).join("");
  const orth = orthologInOther(row);
  const orthBtn = orth ? `<button class="orthjump" id="orthJump">Same gene in ${ORGANISM_META[orth.org].name}: <i>${esc(orth.gene || orth.acc)}</i></button>` : "";
  const weighted = showWeighted() && isNum(row.__c)
    ? `<div class="gcweighted">Weighted score (yours) <b>${row.__c.toFixed(2)}</b> · rank #${compRank(row).toLocaleString()} of ${DATA.rows.length.toLocaleString()}</div>` : "";
  $("drawer").innerHTML = `
    <div class="dhead">
      <button class="close" id="drawerClose" aria-label="Close">×</button>
      <div class="cnav"><button id="cardPrev" aria-label="Previous protein" title="Previous (↑)">↑</button><span id="cardPos"></span><button id="cardNext" aria-label="Next protein" title="Next (↓)">↓</button></div>
      <div class="idrow">${starHTML(row.uniprot_accession)}<span class="gene">${esc(row.gene || row.uniprot_accession)}</span>${fcPill}</div>
      <div class="submeta">
        <button class="accbtn" id="accCopy" title="Copy accession">${esc(row.uniprot_accession)} <span class="cpi">⧉</span></button>
        <span class="org"><span class="odot odot-${state.org}"></span>${ORGANISM_META[state.org].name} ${ORGANISM_META[state.org].strain}</span>
        ${isNum(row.length) ? `<span class="orglen" data-tip="Protein length (amino acids)">${row.length.toLocaleString()} aa</span>` : ""}
      </div>
      ${row.name ? `<div class="pdesc">${esc(row.name)}</div>` : ""}
      ${weighted}
      <div class="links">${links}</div>
      ${orthBtn ? `<div class="orthrow">${orthBtn}</div>` : ""}
    </div>
    <div class="dbody">
      <section class="panel viewerpanel"><div class="pdbslot" id="pdbslot"><div class="pdbmsg">Loading 3D structure…</div></div><p class="pdbcap" id="pdbcap"></p></section>
      <div id="groupPanels"></div>
    </div>`;
  renderGroupPanels(row);
  $("drawerClose").onclick = closeDrawer;
  const pos = filtered.indexOf(row);
  $("cardPos").textContent = pos >= 0 ? `${(pos + 1).toLocaleString()} of ${filtered.length.toLocaleString()}` : "";
  $("cardPrev").disabled = pos <= 0; $("cardNext").disabled = pos < 0 || pos >= filtered.length - 1;
  $("cardPrev").onclick = () => stepDrawer(-1); $("cardNext").onclick = () => stepDrawer(1);
  if (!$("drawer").classList.contains("open")) {
    _lastFocus = document.activeElement;
    $("drawer").classList.add("fresh"); setTimeout(() => $("drawer").classList.remove("fresh"), 900);
  }
  const cardStar = $("drawer").querySelector(".idrow .star");
  if (cardStar) cardStar.onclick = () => toggleStar(row.uniprot_accession);
  const oj = $("orthJump");
  if (oj && orth) oj.onclick = () => {
    document.querySelector(`#orgToggle button[data-org="${orth.org}"]`)?.setAttribute("aria-selected", "true");
    loadOrg(orth.org).then(() => { const r = DATA.rows.find((x) => x.uniprot_accession === orth.acc); if (r) openDrawer(r); });
  };
  const accBtn = $("accCopy");
  if (accBtn) accBtn.onclick = () => {
    try { navigator.clipboard && navigator.clipboard.writeText(row.uniprot_accession); } catch (e) {}
    accBtn.classList.add("copied"); setTimeout(() => accBtn.classList.remove("copied"), 1100);
  };
  $("drawer").classList.add("open");
  $("drawer").setAttribute("aria-hidden", "false");
  $("scrim").classList.add("open");
  renderStructure(row);
}
let _lastFocus = null;
// move through the filtered list without closing the card
function stepDrawer(dir) {
  const i = filtered.indexOf(DATA.rows.find((r) => r.uniprot_accession === state._sel));
  const next = filtered[i + dir];
  if (i < 0 || !next) return;
  openDrawer(next);   // the table behind stays exactly where it is
}
function closeDrawer() {
  $("drawer").classList.remove("open");
  $("drawer").setAttribute("aria-hidden", "true");
  $("scrim").classList.remove("open");
  state._sel = null;
  if (_lastFocus && document.body.contains(_lastFocus)) { try { _lastFocus.focus({ preventScroll: true }); } catch (e) {} }
  document.querySelectorAll("#tbody tr.sel").forEach((tr) => tr.classList.remove("sel"));
}

// ---------- range sliders (both ends) -------------------------------------
const RANGE_PTS = 26;                       // histogram bins
function rangeSpec(spec) {
  const memo = DATA.__rs || (DATA.__rs = {});
  if (memo[spec.key]) return memo[spec.key];
  let { min, max } = spec;
  if (min === undefined) {
    const v = DATA.rows.map((r) => r[spec.key]).filter(isNum);
    min = Math.min(...v); max = Math.max(...v);
  }
  const s = { ...spec, min, max };
  const toT = (x) => spec.scale === "log" ? (Math.log(Math.max(x, 1)) - Math.log(Math.max(min, 1))) / (Math.log(max) - Math.log(Math.max(min, 1))) : (x - min) / (max - min || 1);
  const fromT = (t) => {
    const v = spec.scale === "log" ? Math.exp(Math.log(Math.max(min, 1)) + t * (Math.log(max) - Math.log(Math.max(min, 1)))) : min + t * (max - min);
    const st = spec.step || 0.01;
    return +Math.max(min, Math.min(max, Math.round(v / st) * st)).toFixed(4);
  };
  s.toT = toT; s.fromT = fromT;
  const bins = new Array(RANGE_PTS).fill(0);
  for (const r of DATA.rows) { const x = r[spec.key]; if (isNum(x) && x >= min - 1e-9) bins[Math.min(RANGE_PTS - 1, Math.max(0, Math.floor(toT(x) * RANGE_PTS)))]++; }
  // A bin holding most proteins (ligandability: ~97% pinned at 0, "nobody looked") would flatten every other bar,
  // and no transform can spread a tie. So scale the bars to the rest and draw that bin full height, marked as cut off.
  const total = bins.reduce((a, b) => a + b, 0);
  s.spike = bins.findIndex((c) => c > total / 2);
  const top = Math.max(...bins.filter((_, i) => i !== s.spike).map(Math.sqrt)) || 1;
  s.hist = bins.map((c, i) => i === s.spike ? 1 : Math.sqrt(c) / top);   // square root so tall bins do not hide the rest
  s.counts = bins;
  return (memo[spec.key] = s);
}
const fmtRange = (s, x) => (s.step >= 1 ? Math.round(x).toLocaleString() : x.toFixed(2));
function makeRange(spec0, compact) {
  const s = rangeSpec(spec0);
  const cur = state.filters.ranges[s.key] || [s.min, s.max];
  let lo = cur[0], hi = cur[1];
  const el = document.createElement("div");
  el.className = "rs" + (compact ? " compact" : "");
  el.dataset.key = s.key;
  el.style.setProperty("--c", s.color || "var(--brand)");
  const bars = compact ? "" : `<div class="rs-hist">${s.hist.map((h, i) => i === s.spike
    ? `<i class="spike" style="height:100%;--b:100%" title="${s.counts[i].toLocaleString()} proteins (bar cut off)"></i>`
    : `<i style="height:${Math.max(4, Math.round(h * 100))}%;--b:100%"></i>`).join("")}</div>`;
  el.innerHTML = `<div class="rs-head"><span class="rs-label">${s.label}</span><span class="rs-val"></span></div>${bars}`
    + `<div class="rs-track"><div class="rs-fill"></div><button class="rs-thumb" data-w="lo" role="slider" aria-label="${s.label} minimum"></button><button class="rs-thumb" data-w="hi" role="slider" aria-label="${s.label} maximum"></button></div>`;
  const track = el.querySelector(".rs-track"), fill = el.querySelector(".rs-fill"), val = el.querySelector(".rs-val");
  const tl = el.querySelector('[data-w="lo"]'), th = el.querySelector('[data-w="hi"]');
  const hbars = [...el.querySelectorAll(".rs-hist i")];
  const paint = () => {
    const a = s.toT(lo), b = s.toT(hi);
    tl.style.left = a * 100 + "%"; th.style.left = b * 100 + "%";
    fill.style.left = a * 100 + "%"; fill.style.width = (b - a) * 100 + "%";
    const active = lo > s.min + 1e-9 || hi < s.max - 1e-9;
    el.classList.toggle("active", active);
    val.textContent = active ? `${fmtRange(s, lo)} – ${fmtRange(s, hi)}` : "any";
    hbars.forEach((e, i) => e.classList.toggle("out", active && ((i + 0.5) / RANGE_PTS < a - 0.02 || (i + 0.5) / RANGE_PTS > b + 0.02)));
    tl.setAttribute("aria-valuenow", lo); th.setAttribute("aria-valuenow", hi);
    tl.setAttribute("aria-valuemin", s.min); tl.setAttribute("aria-valuemax", hi);
    th.setAttribute("aria-valuemin", lo); th.setAttribute("aria-valuemax", s.max);
    return active;
  };
  let raf = 0;
  const push = (final) => {
    const active = paint();
    if (active) state.filters.ranges[s.key] = [lo, hi]; else delete state.filters.ranges[s.key];
    if (s.onStore) s.onStore();
    if (final) { cancelAnimationFrame(raf); raf = 0; save(); recompute(); return; }
    if (!raf) raf = requestAnimationFrame(() => { raf = 0; recompute({ quiet: true }); updateAccCounts(); });
  };
  const move = (which, x, final) => {
    const r = track.getBoundingClientRect();
    const v = s.fromT(Math.max(0, Math.min(1, (x - r.left) / r.width)));
    if (which === "lo") lo = Math.min(v, hi); else hi = Math.max(v, lo);
    push(final);
  };
  const drag = (e, which) => {
    e.preventDefault();
    const w = which || (Math.abs(e.clientX - tl.getBoundingClientRect().left) <= Math.abs(e.clientX - th.getBoundingClientRect().left) ? "lo" : "hi");
    const target = w === "lo" ? tl : th;
    target.setPointerCapture?.(e.pointerId); target.focus({ preventScroll: true }); el.classList.add("drag");
    if (!which) move(w, e.clientX);
    const mv = (ev) => move(w, ev.clientX);
    const up = (ev) => { target.removeEventListener("pointermove", mv); target.removeEventListener("pointerup", up); target.removeEventListener("pointercancel", up); el.classList.remove("drag"); move(w, ev.clientX, true); };
    target.addEventListener("pointermove", mv); target.addEventListener("pointerup", up); target.addEventListener("pointercancel", up);
  };
  tl.addEventListener("pointerdown", (e) => drag(e, "lo"));
  th.addEventListener("pointerdown", (e) => drag(e, "hi"));
  track.addEventListener("pointerdown", (e) => { if (!e.target.classList.contains("rs-thumb")) drag(e, null); });
  const key = (which) => (e) => {
    const st = (s.step || 0.01) * (e.shiftKey ? 10 : 1);
    let d = 0;
    if (e.key === "ArrowLeft" || e.key === "ArrowDown") d = -st; else if (e.key === "ArrowRight" || e.key === "ArrowUp") d = st;
    else if (e.key === "Home") { which === "lo" ? (lo = s.min) : (hi = lo); push(true); e.preventDefault(); return; }
    else if (e.key === "End") { which === "hi" ? (hi = s.max) : (lo = hi); push(true); e.preventDefault(); return; }
    else return;
    e.preventDefault();
    if (s.scale === "log") { const t = s.toT(which === "lo" ? lo : hi) + Math.sign(d) * 0.02 * (e.shiftKey ? 5 : 1); const v = s.fromT(Math.max(0, Math.min(1, t))); which === "lo" ? (lo = Math.min(v, hi)) : (hi = Math.max(v, lo)); }
    else if (which === "lo") lo = Math.max(s.min, Math.min(hi, +(lo + d).toFixed(4))); else hi = Math.min(s.max, Math.max(lo, +(hi + d).toFixed(4)));
    push(true);
  };
  tl.addEventListener("keydown", key("lo")); th.addEventListener("keydown", key("hi"));
  el.querySelector(".rs-label").addEventListener("dblclick", () => { lo = s.min; hi = s.max; push(true); });   // double-click the name: clear this filter
  paint();
  el._range = () => [lo, hi];
  return el;
}
// Localization as the cell itself: concentric rings, centre to outside. Click a ring (or a row) to select.
const LOC_ORDER = ["cytoplasm", "cytoplasmic_membrane", "periplasm", "outer_membrane", "extracellular"];
const LOC_RINGS = { cytoplasm: [0, 22], cytoplasmic_membrane: [22, 29], periplasm: [29, 37], outer_membrane: [37, 43], extracellular: [43, 49] };
const LOC_COLOR = { cytoplasm: "#6bbf59", cytoplasmic_membrane: "#fcbf49", periplasm: "#faa08c", outer_membrane: "#e63946", extracellular: "#e91e8c" };
function buildLocCircle() {
  const host = $("bodyLoc"); if (!host) return;
  host.innerHTML = "";
  if (!has_("localization")) return;
  if (!state.filters.cats.localization) state.filters.cats.localization = new Set();
  const sel = state.filters.cats.localization;
  const ring = (k) => {
    const [i, o] = LOC_RINGS[k];
    return i === 0 ? `<circle class="lr disc" data-k="${k}" r="${o}"/>`
      : `<circle class="lr" data-k="${k}" r="${(i + o) / 2}" stroke-width="${o - i}"/>`;
  };
  host.innerHTML = `<div class="locwrap"><svg class="loccircle" viewBox="-50 -50 100 100" role="group" aria-label="Localization">${LOC_ORDER.slice().reverse().map(ring).join("")}`
    + `<g class="mem" fill="none">${[22, 29, 37, 43].map((r) => `<circle r="${r}"/>`).join("")}</g></svg>`
    + `<div class="loclegend">${LOC_ORDER.map((k) => `<button class="locrow" data-k="${k}" aria-pressed="false"><span>${LOC_LABEL[k]}</span></button>`).join("")}</div></div>`;
  const toggle = (k) => { if (sel.has(k)) sel.delete(k); else sel.add(k); paintLoc(); save(); recompute(); };
  host.querySelectorAll(".lr").forEach((c) => c.addEventListener("click", () => toggle(c.dataset.k)));
  host.querySelectorAll(".locrow").forEach((b) => { b.addEventListener("click", () => toggle(b.dataset.k));
    b.addEventListener("mouseenter", () => host.querySelector(`.lr[data-k="${b.dataset.k}"]`).classList.add("hov"));
    b.addEventListener("mouseleave", () => host.querySelector(`.lr[data-k="${b.dataset.k}"]`).classList.remove("hov")); });
  host.querySelectorAll(".lr").forEach((c) => { c.addEventListener("mouseenter", () => host.querySelector(`.locrow[data-k="${c.dataset.k}"]`).classList.add("hov")); c.addEventListener("mouseleave", () => host.querySelector(`.locrow[data-k="${c.dataset.k}"]`).classList.remove("hov")); });
  paintLoc();
}
// counts are what is left once every OTHER filter applies; rings fade with them
function paintLoc() {
  const host = $("bodyLoc"); if (!host || !host.querySelector(".loccircle")) return;
  const sel = state.filters.cats.localization || new Set();
  const counts = {};
  for (const r of DATA.rows) { if (passFilters(r, "localization")) counts[r.localization] = (counts[r.localization] || 0) + 1; }
  const max = Math.max(1, ...LOC_ORDER.map((k) => counts[k] || 0));
  host.querySelectorAll(".lr").forEach((c) => {
    const k = c.dataset.k, on = sel.has(k), n = counts[k] || 0;
    c.style.opacity = (sel.size ? (on ? 1 : 0.5) : 0.55 + 0.45 * Math.sqrt(n / max)).toFixed(2);   // darker = more proteins left
    c.classList.toggle("on", on);
  });
  host.querySelectorAll(".locrow").forEach((b) => {
    const k = b.dataset.k; b.setAttribute("aria-pressed", String(sel.has(k)));
    b.classList.toggle("empty", (counts[k] || 0) === 0 && !sel.has(k));   // nothing left in this compartment
  });
}
// Function as compact chips, grouped as COG groups them. All ON by default (nothing filtered). The stored filter
// is still the set of INCLUDED categories, empty = all, so saved views keep working.
//   click a chip: turn it off / on      double-click: only this one      click a group name: toggle the group
const chipCodes = (k) => COG_CHIP_CODES[k] || [k];
function cogChips() {   // chips with at least one of their codes in this proteome
  const present = new Set(distinct("cog"));
  return COG_GROUPS.flatMap(([, ks]) => ks).filter((k) => chipCodes(k).some((c) => present.has(c)));
}
function cogPresent() { const present = new Set(distinct("cog")); return cogChips().flatMap(chipCodes).filter((c) => present.has(c)); }
function setCog(next) {   // next: Set of codes to show; all (or none) of them collapses to "no filter"
  const all = cogPresent();
  const sel = state.filters.cats.cog;
  sel.clear();
  if (next.size && all.some((c) => !next.has(c))) next.forEach((c) => sel.add(c));
  paintCog(); updateAccCounts(); save(); recompute();
}
function buildCogGrid() {
  const host = $("bodyCog"); if (!host) return;
  host.innerHTML = "";
  if (!has_("cog")) return;
  if (!state.filters.cats.cog) state.filters.cats.cog = new Set();
  const chipSet = new Set(cogChips()), codes = cogPresent();
  const tip = (k) => chipCodes(k).map((c) => c === "none" ? "no COG hit" : c + ": " + ((FC_BY_ID[c] || {}).label || "")).join(" + ");
  host.innerHTML = `<div class="coggrid">` + COG_GROUPS.map(([g, ks]) => {
    const chips = ks.filter((k) => chipSet.has(k)).map((k) => `<button class="cogt" data-k="${k}" aria-pressed="true" title="${tip(k)} (double-click: only this)">${COG_SHORT[k] || k}</button>`).join("");
    return chips ? `<div class="cogg"><button class="cog-h" data-g="${g}" title="Turn the whole group on or off">${g}</button><div class="cogchips">${chips}</div></div>` : "";
  }).join("") + `</div>`;
  const shown = () => { const sel = state.filters.cats.cog; return new Set(sel.size ? sel : codes); };
  const isOn = (on, k) => chipCodes(k).some((c) => on.has(c));
  host.querySelectorAll(".cogt").forEach((b) => {
    const k = b.dataset.k;
    b.addEventListener("click", () => {
      const on = shown();
      if (isOn(on, k)) { const rest = [...on].filter((c) => !chipCodes(k).includes(c)); if (rest.length) chipCodes(k).forEach((c) => on.delete(c)); }
      else chipCodes(k).forEach((c) => on.add(c));   // the last chip on stays on
      setCog(on);
    });
    b.addEventListener("dblclick", () => setCog(new Set(chipCodes(k))));
  });
  host.querySelectorAll(".cog-h").forEach((h) => h.addEventListener("click", () => {
    const ks = (COG_GROUPS.find(([g]) => g === h.dataset.g) || [, []])[1].filter((k) => chipSet.has(k));
    const on = shown(), allOn = ks.every((k) => isOn(on, k));
    ks.forEach((k) => chipCodes(k).forEach((c) => (allOn ? on.delete(c) : on.add(c))));
    setCog(on.size ? on : new Set(ks.flatMap(chipCodes)));
  }));
  paintCog();
}
function paintCog() {
  const host = $("bodyCog"); if (!host || !host.querySelector(".coggrid")) return;
  const sel = state.filters.cats.cog || new Set();
  const counts = {};
  for (const r of DATA.rows) if (passFilters(r, "cog")) counts[r.cog] = (counts[r.cog] || 0) + 1;
  host.querySelectorAll(".cogt").forEach((b) => {
    const ks = chipCodes(b.dataset.k), on = !sel.size || ks.some((c) => sel.has(c));
    b.setAttribute("aria-pressed", String(on));
    b.classList.toggle("empty", ks.every((c) => !counts[c]));
  });
  host.querySelectorAll(".cog-h").forEach((h) => {
    const ks = [...h.nextElementSibling.querySelectorAll(".cogt")];
    h.classList.toggle("off", ks.every((b) => b.getAttribute("aria-pressed") === "false"));
  });
}
// Each histogram shows how many proteins are left in every bin once ALL the other filters are applied (dark),
// against the full distribution (light), so a slider tells you what the rest of the filters have already removed.
function refreshHists() {
  if (!DATA) return;
  for (const el of document.querySelectorAll("#rangePanel .rs")) {
    const s = rangeSpec(RANGE_FILTERS.find((f) => f.key === el.dataset.key) || {});
    if (!s.counts) continue;
    const bins = new Array(RANGE_PTS).fill(0);
    for (const r of DATA.rows) {
      if (!passFilters(r, s.key)) continue;
      const x = r[s.key];
      if (isNum(x) && x >= s.min - 1e-9) bins[Math.min(RANGE_PTS - 1, Math.max(0, Math.floor(s.toT(x) * RANGE_PTS)))]++;
    }
    const bs = el.querySelectorAll(".rs-hist i");
    for (let i = 0; i < bs.length; i++) bs[i].style.setProperty("--b", s.counts[i] ? Math.round(Math.sqrt(bins[i] / s.counts[i]) * 100) + "%" : "0%");
  }
}
function refreshAll() { refreshHists(); paintLoc(); paintCog(); }
const RANGE_SPECS = () => RANGE_FILTERS.filter((f) => DATA.columns.includes(f.key));
function buildRangePanel() {
  keepSidebarScroll(() => {
    const host = $("rangePanel"); if (!host) return;
    host.innerHTML = "";
    for (const f of RANGE_SPECS()) host.appendChild(makeRange(f, false));
    updateAccCounts();
    refreshHists();
  });
}

// ---------- sidebar builders ----------------------------------------------
// ---- weighted score as a radar ---------------------------------------------------------
// One spoke per axis. A weight is the handle's distance from the centre (0-100, steps of 5); the centre is
// "off". Weights are independent: the score divides by their sum, and each tip shows its share in percent.
const WSTEP = 5;
// legacy state: `on:false` meant "ignore this weight"; now a weight of 0 says that
function tidyWeights() {
  for (const c of COMPONENTS) {
    const w = state.weights[c.key] || { on: c.on, weight: c.weight };
    let v = w.on && isNum(w.weight) ? w.weight : 0;
    v = Math.max(0, Math.min(100, Math.round(v / WSTEP) * WSTEP));
    state.weights[c.key] = { on: v > 0, weight: v };
  }
}
function setWeights(map) {   // { key: 0-100 }; missing keys -> 0
  for (const c of COMPONENTS) { const v = map[c.key] || 0; state.weights[c.key] = { on: v > 0, weight: v }; }
  weightsDirty = true; buildWeights(); save(); recompute();
}
const DEFAULT_WEIGHTS = () => Object.fromEntries(COMPONENTS.map((c) => [c.key, c.weight]));
function buildWeights() {
  const host = $("weights"); if (!host) return;
  host.innerHTML = "";
  tidyWeights();
  const axes = COMPONENTS.filter((c) => AVAIL[c.key]);
  if (!axes.length) return;
  const W = 260, H = 226, cx = 130, cy = 112, R = 70, n = axes.length;
  const ang = (i) => -Math.PI / 2 + i * 2 * Math.PI / n;
  const at = (i, r) => [cx + r * Math.cos(ang(i)), cy + r * Math.sin(ang(i))];
  const f = (x) => x.toFixed(1);
  const rings = [0.25, 0.5, 0.75, 1].map((t) => `<circle class="wr-ring" cx="${cx}" cy="${cy}" r="${R * t}"/>`).join("");
  const spokes = axes.map((c, i) => { const [x, y] = at(i, R); return `<line class="wr-spoke" x1="${cx}" y1="${cy}" x2="${f(x)}" y2="${f(y)}"/>`; }).join("");
  const labels = axes.map((c, i) => {
    const [x, y] = at(i, R + 16), cos = Math.cos(ang(i)), sin = Math.sin(ang(i));
    const anchor = cos > 0.3 ? "start" : cos < -0.3 ? "end" : "middle";
    const dy = sin < -0.5 ? -14 : sin > 0.5 ? 6 : -4;
    const name = `<tspan x="${f(x)}">${c.radar || c.label}</tspan>`;   // one line: long names carry a short `radar` form
    return `<text class="wr-lab" data-i="${i}" x="${f(x)}" y="${f(y + dy)}" text-anchor="${anchor}"${tipAttr(c.help + " Click to switch on or off.")}>${name}`
      + `<tspan class="wr-pct" x="${f(x)}" dy="13" style="fill:${colColor(c.key)}"></tspan></text>`;
  }).join("");
  const handles = axes.map((c, i) => `<circle class="wr-h" data-i="${i}" r="6.5" tabindex="0" role="slider" aria-label="${c.label} weight" aria-valuemin="0" aria-valuemax="100" style="--hc:${colColor(c.key)}"/>`).join("");
  host.innerHTML = `<svg class="wradar" viewBox="0 0 ${W} ${H}" role="group" aria-label="Weighted score: drag each dot along its spoke">`
    + rings + spokes + `<polygon class="wr-poly"/>` + `<text class="wr-empty" x="${cx}" y="${cy + 4}" text-anchor="middle">drag a spoke</text>`
    + handles + labels + `</svg>`;
  const svg = host.querySelector("svg"), poly = svg.querySelector(".wr-poly"), empty = svg.querySelector(".wr-empty");
  const hs = [...svg.querySelectorAll(".wr-h")], pcts = [...svg.querySelectorAll(".wr-pct")];
  const wOf = (i) => state.weights[axes[i].key].weight;
  const paint = () => {
    const tot = axes.reduce((s, c, i) => s + wOf(i), 0);
    poly.setAttribute("points", axes.map((c, i) => at(i, R * wOf(i) / 100).map(f).join(",")).join(" "));
    poly.style.display = tot ? "" : "none";
    empty.style.display = tot ? "none" : "";
    hs.forEach((h, i) => {
      const v = wOf(i), [x, y] = at(i, v ? R * v / 100 : 9);
      h.setAttribute("cx", f(x)); h.setAttribute("cy", f(y));
      h.classList.toggle("zero", !v);
      h.setAttribute("aria-valuenow", v);
      pcts[i].textContent = tot && v ? Math.round(100 * v / tot) + "%" : "off";
      pcts[i].classList.toggle("off", !v);
    });
  };
  const set = (i, v) => { v = Math.max(0, Math.min(100, Math.round(v / WSTEP) * WSTEP)); state.weights[axes[i].key] = { on: v > 0, weight: v }; paint(); };
  const commit = () => { weightsDirty = true; save(); recompute(); };
  // pointer -> (spoke, weight): the spoke is the handle grabbed, else the nearest by angle
  const local = (e) => { const p = new DOMPoint(e.clientX, e.clientY).matrixTransform(svg.getScreenCTM().inverse()); return [p.x - cx, p.y - cy]; };
  const nearest = (dx, dy) => { let best = 0, bd = Infinity; axes.forEach((c, i) => { const d = Math.abs(Math.atan2(Math.sin(Math.atan2(dy, dx) - ang(i)), Math.cos(Math.atan2(dy, dx) - ang(i)))); if (d < bd) { bd = d; best = i; } }); return best; };
  const valueOn = (i, dx, dy) => 100 * Math.max(0, dx * Math.cos(ang(i)) + dy * Math.sin(ang(i))) / R;
  let drag = null;
  svg.addEventListener("pointerdown", (e) => {
    if (e.target.closest(".wr-lab")) return;
    const [dx, dy] = local(e);
    if (Math.hypot(dx, dy) > R + 14) return;
    const h = e.target.closest(".wr-h");
    drag = h ? +h.dataset.i : nearest(dx, dy);
    svg.setPointerCapture(e.pointerId); svg.classList.add("dragging");
    set(drag, valueOn(drag, dx, dy)); e.preventDefault();
  });
  svg.addEventListener("pointermove", (e) => { if (drag === null) return; const [dx, dy] = local(e); set(drag, valueOn(drag, dx, dy)); });
  const end = () => { if (drag === null) return; drag = null; svg.classList.remove("dragging"); commit(); };
  svg.addEventListener("pointerup", end); svg.addEventListener("pointercancel", end);
  hs.forEach((h, i) => h.addEventListener("keydown", (e) => {
    const d = { ArrowUp: 1, ArrowRight: 1, ArrowDown: -1, ArrowLeft: -1 }[e.key];
    if (!d) return; e.preventDefault(); set(i, wOf(i) + d * WSTEP); commit();
  }));
  svg.querySelectorAll(".wr-lab").forEach((t) => t.addEventListener("click", () => {
    const i = +t.dataset.i; set(i, wOf(i) ? 0 : (axes[i].weight || 50)); commit();
  }));
  paint();
}
function distinct(key) {
  const s = new Set();
  for (const r of DATA.rows) if (r[key] !== null && r[key] !== undefined && r[key] !== "") s.add(r[key]);
  return [...s].sort();
}
// "No human ortholog" as one switch: on keeps only proteins without one (a human ortholog is a
// selectivity liability), off shows everything. A saved "yes" from the old any/yes/no control reads as off.
function noHumanSwitch() {
  const key = "has_human_ortholog";
  if (state.filters.bools[key] !== "no") delete state.filters.bools[key];
  const btn = document.createElement("button");
  btn.className = "switch";
  btn.setAttribute("role", "switch");
  btn.setAttribute("aria-label", "Only proteins with no human ortholog");
  btn.title = "Hide proteins with a human ortholog (OrthoFinder or reciprocal best hit)";
  const sync = () => btn.setAttribute("aria-checked", String(state.filters.bools[key] === "no"));
  const flip = () => {
    if (state.filters.bools[key] === "no") delete state.filters.bools[key]; else state.filters.bools[key] = "no";
    sync(); save(); recompute(); updateAccCounts();
  };
  btn.addEventListener("click", flip);
  const title = document.querySelector("#accOrth .fs-title");
  if (title && !title._wired) { title._wired = true; title.addEventListener("click", () => btn.click()); }
  sync();
  return btn;
}
function anyFilterActive() {
  const f = state.filters;
  return Object.keys(f.ranges).length > 0
    || Object.values(f.cats).some((s) => s && s.size)
    || Object.values(f.bools).some((m) => m && m !== "any");
}
function chipRow(host, key, { label, short, order, text, inline } = {}) {
  if (!state.filters.cats[key]) state.filters.cats[key] = new Set();
  const active = state.filters.cats[key];
  const f = CATEGORICAL_FILTERS.find((x) => x.key === key) || {};
  const row = document.createElement("div");
  row.className = "chiprowf" + (inline ? " inline" : "");
  if (label) row.innerHTML = `<span class="flabel">${label}</span>`;
  const chk = document.createElement("div"); chk.className = "chk"; row.appendChild(chk);
  const present = new Set(distinct(key));
  const values = (order || f.order) ? (order || f.order).filter((v) => present.has(v)) : [...present];
  for (const val of values) {
    const meta = f.meta && f.meta[val];
    const b = document.createElement("button");
    b.textContent = text ? text(val) : (short && meta ? meta.abbr : (meta ? meta.label : String(val).replace(/_/g, " ")));
    if (meta) b.title = meta.label;
    b.setAttribute("aria-pressed", active.has(val) ? "true" : "false");
    if (meta) { b.style.setProperty("--cc", meta.color); b.classList.add("cchip"); }
    b.onclick = () => {
      if (active.has(val)) active.delete(val); else active.add(val);
      b.setAttribute("aria-pressed", active.has(val) ? "true" : "false");
      save(); recompute();
    };
    chk.appendChild(b);
  }
  host.appendChild(row);
}
function keepSidebarScroll(fn) {
  const sb = $("sidebar"), y = sb ? sb.scrollTop : 0;
  fn();
  if (sb) sb.scrollTop = y;
}
function buildCatPanel() { keepSidebarScroll(buildCatPanelNow); }
function buildCatPanelNow() {
  const fill = (id) => { const e = $(id); if (e) e.innerHTML = ""; return e; };
  buildLocCircle();
  const og = fill("bodyOrth"); if (og && has_("has_human_ortholog")) og.appendChild(noHumanSwitch());
  buildCogGrid();
  updateAccCounts();
}
// the little count on a section = how many selections are active inside it
function updateAccCounts() {
  const set = (id, n) => { const e = document.querySelector(`#${id} .acount`); if (e) e.textContent = n ? String(n) : ""; };
  const size = (k) => (state.filters.cats[k] ? state.filters.cats[k].size : 0);
  set("accLoc", size("localization"));
  const cogSel = state.filters.cats.cog, cogT = document.querySelector("#accCog .acount");
  const cogN = cogSel && cogSel.size ? cogChips().filter((k) => chipCodes(k).some((c) => cogSel.has(c))).length : 0;
  if (cogT) cogT.innerHTML = cogN ? `${cogN} of ${cogChips().length} <button class="cogall" type="button">show all</button>` : "";
  const ca = cogT && cogT.querySelector(".cogall"); if (ca) ca.onclick = () => setCog(new Set());
  set("accOrth", (state.filters.bools.has_human_ortholog || "any") !== "any" ? 1 : 0);
  const b = $("sbClear"); if (b) b.hidden = !(anyFilterActive() || !!state.filters.search);
}
function buildFilterBar() { buildCatPanel(); }
// ---------- presets --------------------------------------------------------
// "reset" is the only preset left: clear every filter
function applyPreset() {
  state.filters.search = ""; if ($("search")) $("search").value = "";
  state.filters.ranges = {};
  state.filters.cats = {};
  state.filters.bools = {};
  state.mapSel = null;
  buildRangePanel();
  buildFilterBar();
  save(); recompute();
}
// same-gene ortholog in the other organism (gene-symbol match; needs the other org cached)
function orthologInOther(row) {
  const other = state.org === "kp" ? "ec" : "kp";
  const data = cache[other];
  const g = (row.gene || "").trim().toLowerCase();   // same gene symbol only; not an orthology call
  if (!data || !g) return null;
  const m = data.rows.find((r) => (r.gene || "").trim().toLowerCase() === g);
  return m ? { org: other, acc: m.uniprot_accession, gene: m.gene } : null;
}

// ---------- CSV export -----------------------------------------------------
function exportCSV() {
  const cols = visibleColumns();
  const header = ["rank", "uniprot_accession", "gene", ...cols.map((c) => c.key)];   // the weighted score is never exported
  const lines = [header.join(",")];
  filtered.forEach((row, i) => {
    const vals = [i + 1, row.uniprot_accession, row.gene || ""];
    for (const c of cols) { let v = row[c.key]; v = (v === null || v === undefined) ? "" : v; vals.push(v); }
    lines.push(vals.map((v) => {
      const s = String(v);
      return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
    }).join(","));
  });
  const blob = new Blob([lines.join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `gradi_${state.org === "kp" ? "kpneumoniae" : "ecoli"}_targets.csv`;
  a.click();
  URL.revokeObjectURL(a.href);
}

// ---------- organism loading ----------------------------------------------
// Convert the columnar wire format (rows = arrays of cell values in `columns`
// order) back into the array-of-objects shape the rest of the app assumes.
// Falls back to the payload unchanged if it's already row-objects.
function rehydrate(payload) {
  if (!payload || payload.format !== "columnar") return payload;
  const cols = payload.columns, n = cols.length;
  payload.rows = payload.rows.map((arr) => {
    const o = {};
    for (let j = 0; j < n; j++) o[cols[j]] = arr[j];
    return o;
  });
  delete payload.format;
  return payload;
}
// Warm the cache for the other organism during idle time so the toggle is instant.
function prefetchOther(org) {
  const other = org === "kp" ? "ec" : "kp";
  if (cache[other] || !ORGANISM_META[other]) return;
  const run = () => {
    if (cache[other]) return;
    fetch(ORGANISM_META[other].file)
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => { if (j && !cache[other]) cache[other] = rehydrate(j); })
      .catch(() => {});
  };
  if (window.requestIdleCallback) requestIdleCallback(run, { timeout: 3000 });
  else setTimeout(run, 1200);
}
const SKELETON = Array.from({ length: 9 }, () => `<div class="skel"><i></i><i></i><i></i><i></i></div>`).join("");
async function loadOrg(org) {
  state.org = org;
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.setAttribute("aria-selected", b.dataset.org === org));
  syncSeg($("orgToggle"));
  if (!cache[org]) { $("table").hidden = true; $("loading").hidden = false; $("loading").innerHTML = SKELETON; }
  if (!cache[org]) {
    const meta = ORGANISM_META[org];
    const res = await fetch(meta.file);
    if (!res.ok) { $("loading").textContent = `Failed to load ${meta.file} (${res.status})`; return; }
    cache[org] = rehydrate(await res.json());
  }
  DATA = cache[org];
  weightsDirty = true;   // composite must be (re)computed for the newly-active organism
  _pctCache = {};        // percentile distributions are per-organism
  AVAIL = {};
  for (const c of COMPONENTS) AVAIL[c.key] = DATA.columns.includes(c.key) && (c.available !== false);
  // seed visible columns from the active table view (fall back to overview)
  $("genInfo").textContent = `${ORGANISM_META[org].name} ${ORGANISM_META[org].strain} · UniProt ${DATA.proteome_id} · exported ${DATA.generated_at}`;
  buildViewSwitch();
  buildWeights(); buildRangePanel();
  buildFilterBar();
  stashEssentiality(); applyTransferMode();
  closeDrawer();
  recompute();
  save();
  prefetchOther(org);
}

// ---------- hover tooltips ([data-tip]) ------------------------------------
let _tipEl = null;
function initTooltips() {
  _tipEl = document.createElement("div");
  _tipEl.className = "hovertip";
  _tipEl.hidden = true;
  document.body.appendChild(_tipEl);
  const show = (target, x, y) => {
    _tipEl.textContent = target.getAttribute("data-tip");
    _tipEl.hidden = false;
    const w = _tipEl.offsetWidth, vw = innerWidth;
    let left = x + 14; if (left + w + 10 > vw) left = x - w - 14;
    _tipEl.style.left = Math.max(8, left) + "px";
    _tipEl.style.top = (y + 16) + "px";
  };
  document.addEventListener("mousemove", (e) => {
    const t = e.target.closest("[data-tip]");
    if (t) show(t, e.clientX, e.clientY);
    else _tipEl.hidden = true;
  });
  document.addEventListener("mouseleave", () => { _tipEl.hidden = true; });
}

// ---------- projection map -------------------------------------------------
let mapPts = [];          // {sx, sy, row} for passing points (hit-testing)
let _hoverAcc = null;
let _mapDrag = null;      // active box-select rectangle {x0,y0,x1,y1}
let _mapDragged = false;  // did the last interaction move enough to be a drag?
function cssColor(varRef) {
  const p = document.createElement("span");
  p.style.cssText = "position:absolute;visibility:hidden;color:" + varRef;
  document.body.appendChild(p);
  const m = getComputedStyle(p).color.match(/(\d+(\.\d+)?)/g);
  p.remove();
  return m ? { r: +m[0], g: +m[1], b: +m[2] } : { r: 120, g: 120, b: 120 };
}
// Default "Show >=" for a colour-by field, from the data so it suits each field and each species:
// ranks are uniform by construction, so their top fifth (the 80th percentile); a field that is zero for most
// proteins (ligands) shows everything above zero; a field piled up at its maximum (cytoplasmic fraction) shows the maximum.
function defaultMapThr(key) {
  const memo = DATA.__mt || (DATA.__mt = {});
  if (memo[key] !== undefined) return memo[key];
  const v = DATA.rows.map((r) => r[key]).filter(isNum).sort((x, y) => x - y);
  let t = 0;
  if (v.length) {
    const zeros = v.filter((x) => x <= 0).length / v.length;
    if (zeros > 0.5) t = 0.01;
    else t = Math.round(v[Math.min(v.length - 1, Math.floor(0.8 * v.length))] * 100) / 100;
  }
  return (memo[key] = Math.max(0, Math.min(1, t)));
}
function mapThr(key) { return state.mapThr[key] !== undefined ? state.mapThr[key] : defaultMapThr(key); }
function syncMapSlider(catMode) {
  const thr = $("mapThresh"); if (!thr) return;
  const key = state.mapColorBy;
  const cat = catMode !== undefined ? catMode : (key === "family" || !!MAP_CATEGORICAL[key]);
  thr.disabled = cat;
  thr.value = cat ? 0 : mapThr(key);
  $("mapThreshVal").textContent = (+thr.value).toFixed(2);
  thr.title = cat ? "" : "Double-click to reset to the default for this view";
  const lt = document.querySelector(".maplabtog"); if (lt) lt.style.visibility = cat ? "visible" : "hidden";   // labels only apply to categories
}
// the weighted score is only offered as a colouring while it exists
let _mapSelSig = null;
function buildMapSelect() {
  const sel = $("mapColor"); if (!sel) return;
  const sig = showWeighted() ? "w" : "n";
  if (sig === _mapSelSig) return;
  _mapSelSig = sig;
  if (state.mapColorBy === "__c" && sig === "n") state.mapColorBy = "degradability_consensus";
  sel.innerHTML = MAP_COLORS.filter((m) => m.key !== "__c" || sig === "w").map((m) => `<option value="${m.key}">${m.label}</option>`).join("");
  sel.value = state.mapColorBy;
}
function mapExtent() {
  if (DATA.__extent) return DATA.__extent;
  let minx = Infinity, maxx = -Infinity, miny = Infinity, maxy = -Infinity;
  for (const r of DATA.rows) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    if (r.tsne_x < minx) minx = r.tsne_x; if (r.tsne_x > maxx) maxx = r.tsne_x;
    if (r.tsne_y < miny) miny = r.tsne_y; if (r.tsne_y > maxy) maxy = r.tsne_y;
  }
  return (DATA.__extent = { minx, maxx, miny, maxy });
}
function drawMap() {
  buildMapSelect();
  const wrap = $("mapwrap"), base = $("mapcanvas"), over = $("mapoverlay");
  const W = wrap.clientWidth, H = wrap.clientHeight, dpr = devicePixelRatio || 1;
  for (const c of [base, over]) { c.width = W * dpr; c.height = H * dpr; }
  const ctx = base.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);

  const { minx, maxx, miny, maxy } = mapExtent();
  const P = 26, sx = (W - 2 * P) / ((maxx - minx) || 1), sy = (H - 2 * P) / ((maxy - miny) || 1);
  const s = Math.min(sx, sy);
  const ox = P + (W - 2 * P - s * (maxx - minx)) / 2, oy = P + (H - 2 * P - s * (maxy - miny)) / 2;
  const X = (x) => ox + (x - minx) * s;
  const Y = (y) => oy + (maxy - y) * s;   // flip

  const colorBy = state.mapColorBy;
  const acc = cssColor(colColor(colorBy)), faint = cssColor("var(--faint)");
  const catMap = MAP_CATEGORICAL[colorBy];                 // named-categorical (or undefined)
  const catMode = colorBy === "family" || !!catMap;
  const L = MAP_CAT_PALETTE.length;
  const catRGB = {};                                       // resolve category CSS vars -> rgb once
  if (catMap) for (const k in catMap) catRGB[k] = cssColor(catMap[k]);
  const t = catMode ? 0 : mapThr(colorBy);
  const rad = (r) => 2 + 3 * Math.max(0, Math.min(1, isNum(r.__c) ? r.__c : 0));   // size ∝ composite
  syncMapSlider(catMode);

  // base layer: every protein as faint context; also collect hit-test points
  mapPts = [];
  ctx.fillStyle = `rgba(${faint.r},${faint.g},${faint.b},0.16)`;
  for (const r of DATA.rows) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    const px = X(r.tsne_x), py = Y(r.tsne_y);
    ctx.beginPath(); ctx.arc(px, py, 1.3, 0, 6.283); ctx.fill();
    mapPts.push({ sx: px, sy: py, row: r });
  }
  // highlight layer: passing rows, coloured by the chosen field
  let hi = 0;
  ctx.lineWidth = 0.6; ctx.strokeStyle = "rgba(255,255,255,0.85)";
  for (const r of filtered) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    let fill;
    if (colorBy === "family") {
      if (!isNum(r.family)) continue;
      fill = MAP_CAT_PALETTE[((r.family % L) + L) % L];
    } else if (catMap) {
      const c = catRGB[r[colorBy]];
      if (!c) continue;
      fill = `rgb(${c.r},${c.g},${c.b})`;
    } else {
      const v = r[colorBy];
      if (!isNum(v) || v < t) continue;
      const vv = Math.max(0, Math.min(1, v));
      fill = `rgba(${acc.r},${acc.g},${acc.b},${0.5 + 0.5 * vv})`;
    }
    ctx.fillStyle = fill;
    ctx.beginPath(); ctx.arc(X(r.tsne_x), Y(r.tsne_y), rad(r), 0, 6.283); ctx.fill(); ctx.stroke();
    hi++;
  }
  // cluster labels: name each category at its median position (categorical + toggle on)
  if (catMap && state.mapLabels !== false) {
    const groups = {};
    for (const r of filtered) {
      const v = r[colorBy];
      if (!isNum(r.tsne_x) || !isNum(r.tsne_y) || v == null || v === "") continue;
      (groups[v] || (groups[v] = [])).push([X(r.tsne_x), Y(r.tsne_y)]);
    }
    const med = (a) => { a.sort((p, q) => p - q); return a[a.length >> 1]; };
    ctx.save();
    ctx.font = '600 11.5px "Inter",system-ui,sans-serif';
    ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.lineJoin = "round";
    for (const v in groups) {
      const pts = groups[v];
      if (pts.length < 12) continue;
      const mx = med(pts.map((p) => p[0])), my = med(pts.map((p) => p[1]));
      const label = catLabel(colorBy, v);
      ctx.lineWidth = 3.5; ctx.strokeStyle = "rgba(255,255,255,0.92)"; ctx.strokeText(label, mx, my);
      ctx.fillStyle = "#2C3E50"; ctx.fillText(label, mx, my);
    }
    ctx.restore();
  }
  const cnt = $("mapHiCount"); if (cnt) cnt.textContent = `${hi.toLocaleString()} ${catMode ? "shown" : "highlighted"}`;
  drawMapOverlay();
  renderMapLegend();
}
function drawMapOverlay() {
  const over = $("mapoverlay"), dpr = devicePixelRatio || 1;
  const ctx = over.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, over.width, over.height);
  const acc = cssColor(colColor(state.mapColorBy));
  const mark = (accn, ring) => {
    const p = mapPts.find((q) => q.row.uniprot_accession === accn);
    if (!p) return;
    ctx.beginPath(); ctx.arc(p.sx, p.sy, ring, 0, 6.283);
    ctx.lineWidth = 2; ctx.strokeStyle = `rgb(${acc.r},${acc.g},${acc.b})`; ctx.stroke();
  };
  // ring shortlisted (starred) proteins in amber
  if (state.shortlist && state.shortlist.length) {
    ctx.lineWidth = 2; ctx.strokeStyle = "#FCBF49";
    for (const a of state.shortlist) {
      const p = mapPts.find((q) => q.row.uniprot_accession === a);
      if (p) { ctx.beginPath(); ctx.arc(p.sx, p.sy, 6, 0, 6.283); ctx.stroke(); }
    }
  }
  if (state._sel) mark(state._sel, 7);
  if (_hoverAcc && _hoverAcc !== state._sel) mark(_hoverAcc, 6);
  // gene label next to the ringed selected/hovered point
  const labelPt = (accn) => {
    const p = mapPts.find((q) => q.row.uniprot_accession === accn);
    if (!p) return;
    const g = p.row.gene || p.row.uniprot_accession;
    ctx.save();
    ctx.font = '600 12px "Inter",system-ui,sans-serif';
    ctx.textBaseline = "middle"; ctx.lineJoin = "round";
    ctx.lineWidth = 3.5; ctx.strokeStyle = "rgba(255,255,255,0.92)"; ctx.strokeText(g, p.sx + 10, p.sy);
    ctx.fillStyle = "#2C3E50"; ctx.fillText(g, p.sx + 10, p.sy);
    ctx.restore();
  };
  if (state._sel) labelPt(state._sel);
  else if (_hoverAcc) labelPt(_hoverAcc);
  if (_mapDrag) {
    const { x0, y0, x1, y1 } = _mapDrag;
    ctx.fillStyle = "rgba(108,92,231,0.10)"; ctx.strokeStyle = "rgba(108,92,231,0.9)"; ctx.lineWidth = 1;
    const rx = Math.min(x0, x1), ry = Math.min(y0, y1), rw = Math.abs(x1 - x0), rh = Math.abs(y1 - y0);
    ctx.fillRect(rx, ry, rw, rh); ctx.strokeRect(rx, ry, rw, rh);
  }
}
function renderMapLegend() {
  const lbl = (MAP_COLORS.find((m) => m.key === state.mapColorBy) || {}).label || "Value";
  if (state.mapColorBy === "family") {
    const sw = MAP_CAT_PALETTE.slice(0, 8).map((c) => `<span style="background:${c}"></span>`).join("");
    $("maplegend").innerHTML = `<div class="lt">${lbl}</div><div class="catramp">${sw}</div>`
      + `<div class="sc"><span>ESM-C cluster</span></div>`;
    return;
  }
  const field = state.mapColorBy, catMap = MAP_CATEGORICAL[field];
  if (catMap) {
    const present = new Set();
    for (const r of DATA.rows) { const v = r[field]; if (v != null && v !== "") present.add(v); }
    const active = state.filters.cats[field];
    const items = Object.keys(catMap).filter((k) => present.has(k));
    const rows = items.map((k) => {
      const on = active && active.size ? active.has(k) : true;
      return `<button class="catrow${on ? "" : " off"}" data-cat="${k}">`
        + `<span class="sw" style="background:${catMap[k]}"></span>${catLabel(field, k)}</button>`;
    }).join("");
    $("maplegend").innerHTML = `<div class="lt">${lbl}</div><div class="catlist">${rows}</div>`
      + `<div class="sc"><span>click to isolate a group</span></div>`;
    $("maplegend").querySelectorAll("button[data-cat]").forEach((b) => {
      b.onclick = () => {
        if (!state.filters.cats[field]) state.filters.cats[field] = new Set();
        const s = state.filters.cats[field];
        if (s.has(b.dataset.cat)) s.delete(b.dataset.cat); else s.add(b.dataset.cat);
        if (!s.size) delete state.filters.cats[field];
        save(); recompute();
      };
    });
    return;
  }
  const col = colColor(state.mapColorBy);
  const prov = isProvisional(state.mapColorBy) ? ` <span class="provtag">provisional</span>` : "";
  $("maplegend").innerHTML = `<div class="lt">${lbl}${prov}</div>`
    + `<div class="ramp" style="background:linear-gradient(90deg,color-mix(in srgb,${col} 14%,#fff),${col})"></div>`
    + `<div class="sc"><span>0</span><span>1</span></div>`;
}
function nearestPt(mx, my) {
  let best = null, bd = 100;
  for (const p of mapPts) {
    const d = (p.sx - mx) ** 2 + (p.sy - my) ** 2;
    if (d < bd) { bd = d; best = p; }
  }
  return best;
}
function initMap() {
  // populate colour-by select
  const sel = $("mapColor");
  buildMapSelect();
  sel.onchange = () => { state.mapColorBy = sel.value; save(); if (state.view === "map") drawMap(); else syncMapSlider(); };

  const thr = $("mapThresh");
  thr.oninput = () => { $("mapThreshVal").textContent = (+thr.value).toFixed(2); };
  thr.onchange = () => {
    state.mapThr[state.mapColorBy] = +thr.value;
    $("mapThreshVal").textContent = (+thr.value).toFixed(2);
    save(); if (state.view === "map") drawMap();
  };
  thr.ondblclick = () => { delete state.mapThr[state.mapColorBy]; syncMapSlider(); save(); if (state.view === "map") drawMap(); };   // double-click: back to the default

  const ml = $("mapLabels");
  if (ml) {
    ml.checked = state.mapLabels !== false;
    ml.onchange = () => { state.mapLabels = ml.checked; save(); if (state.view === "map") drawMap(); };
  }

  const wrap = $("mapwrap");
  wrap.addEventListener("mousedown", (e) => {
    if (e.button !== 0) return;
    const rect = wrap.getBoundingClientRect();
    const x = e.clientX - rect.left, y = e.clientY - rect.top;
    _mapDrag = { x0: x, y0: y, x1: x, y1: y }; _mapDragged = false;
  });
  wrap.addEventListener("mouseup", () => {
    if (!_mapDrag) return;
    const d = _mapDrag; _mapDrag = null;
    if (!_mapDragged) { drawMapOverlay(); return; }   // a click, not a drag
    const xmin = Math.min(d.x0, d.x1), xmax = Math.max(d.x0, d.x1);
    const ymin = Math.min(d.y0, d.y1), ymax = Math.max(d.y0, d.y1);
    const sel = new Set();
    for (const p of mapPts) if (p.sx >= xmin && p.sx <= xmax && p.sy >= ymin && p.sy <= ymax) sel.add(p.row.uniprot_accession);
    state.mapSel = sel.size ? sel : null;
    updateMapSelChip();
    recompute();
  });
  wrap.addEventListener("mousemove", (e) => {
    const rect = wrap.getBoundingClientRect();
    if (_mapDrag) {
      _mapDrag.x1 = e.clientX - rect.left; _mapDrag.y1 = e.clientY - rect.top;
      if (Math.hypot(_mapDrag.x1 - _mapDrag.x0, _mapDrag.y1 - _mapDrag.y0) > 4) _mapDragged = true;
      $("maptip").hidden = true; wrap.style.cursor = "crosshair"; drawMapOverlay(); return;
    }
    const p = nearestPt(e.clientX - rect.left, e.clientY - rect.top);
    const acc = p ? p.row.uniprot_accession : null;
    if (acc !== _hoverAcc) { _hoverAcc = acc; drawMapOverlay(); }
    const tip = $("maptip");
    if (p) {
      const r = p.row, v = r[state.mapColorBy];
      const cl = (MAP_COLORS.find((m) => m.key === state.mapColorBy) || {}).label;
      tip.innerHTML = `<div class="g">${r.gene || r.uniprot_accession}</div>`
        + `<div class="a">${r.uniprot_accession}</div>`
        + `<div class="m">${cl}: <b>${isNum(v) ? v.toFixed(2) : (v == null ? "–" : catLabel(state.mapColorBy, v))}</b></div>`;
      tip.hidden = false;
      let lx = p.sx + 14; if (lx + tip.offsetWidth + 12 > wrap.clientWidth) lx = p.sx - tip.offsetWidth - 14;
      tip.style.left = Math.max(6, lx) + "px";
      tip.style.top = Math.max(6, p.sy - 10) + "px";
      wrap.style.cursor = "pointer";
    } else { tip.hidden = true; wrap.style.cursor = "default"; }
  });
  wrap.addEventListener("mouseleave", () => { _hoverAcc = null; $("maptip").hidden = true; drawMapOverlay(); });
  wrap.addEventListener("click", (e) => {
    if (_mapDragged) { _mapDragged = false; return; }   // ignore the click that ends a drag
    const rect = wrap.getBoundingClientRect();
    const p = nearestPt(e.clientX - rect.left, e.clientY - rect.top);
    if (p) openDrawer(p.row);
  });
  let rAF = null;
  addEventListener("resize", () => {
    if (state.view !== "map") return;
    if (rAF) cancelAnimationFrame(rAF);
    rAF = requestAnimationFrame(drawMap);
  });
}

// ---------- init -----------------------------------------------------------
function initCollapsers() {
  const collapsed = (() => { try { return JSON.parse(localStorage.getItem("gradi-tpb-collapsed") || "{}"); } catch (e) { return {}; } })();
  document.querySelectorAll(".section.collapsible").forEach((sec) => {
    const key = sec.dataset.sec;
    if (collapsed[key]) sec.classList.add("collapsed");
    sec.querySelector("h2").addEventListener("click", (e) => {
      if (e.target.closest("button")) return;   // don't toggle when clicking reset
      sec.classList.toggle("collapsed");
      collapsed[key] = sec.classList.contains("collapsed");
      try { localStorage.setItem("gradi-tpb-collapsed", JSON.stringify(collapsed)); } catch (err) {}
    });
  });
}
function init() {
  load();
  applyHash();          // a shared link overrides saved state
  initTooltips();
  initMap();
  initCollapsers();
  (() => {
    let saved = {}; try { saved = JSON.parse(localStorage.getItem("gradi-web2-acc") || "{}"); } catch (e) {}
    document.querySelectorAll("details.facc").forEach((d) => {
      if (saved[d.id] !== undefined) d.open = !!saved[d.id];
      d.addEventListener("toggle", () => { saved[d.id] = d.open; try { localStorage.setItem("gradi-web2-acc", JSON.stringify(saved)); } catch (e) {} });
    });
  })();
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.onclick = () => { loadOrg(b.dataset.org); });
  $("search").value = state.filters.search;
  let searchTimer = null;
  $("search").oninput = (e) => {
    state.filters.search = e.target.value.trim();
    if (searchTimer) clearTimeout(searchTimer);
    searchTimer = setTimeout(() => { save(); recompute(); }, 70);
  };
  $("resetWeights").onclick = () => setWeights(DEFAULT_WEIGHTS());
  $("sbClear").onclick = () => applyPreset("reset");
  $("exportBtn").onclick = exportCSV;
  $("filtersBtn").onclick = () => document.body.classList.toggle("sbopen");
  $("sbScrim").onclick = () => document.body.classList.remove("sbopen");
  $("sbDone").onclick = () => document.body.classList.remove("sbopen");
  $("methodsBtn").onclick = openMethods;
  $("methodsScrim").onclick = closeMethods;
  updateShortlistBtn();
  $("shortlistBtn").onclick = () => { state.shortlistOnly = !state.shortlistOnly; updateShortlistBtn(); save(); recompute(); };
  $("scrim").onclick = closeDrawer;
  let _tx = null, _ty = null;
  $("drawer").addEventListener("touchstart", (e) => { _tx = e.touches[0].clientX; _ty = e.touches[0].clientY; }, { passive: true });
  $("drawer").addEventListener("touchend", (e) => {
    if (_tx === null) return;
    const dx = e.changedTouches[0].clientX - _tx, dy = e.changedTouches[0].clientY - _ty;
    _tx = null;
    if (dx > 90 && Math.abs(dy) < 50) closeDrawer();
  }, { passive: true });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") { closeDrawer(); closeMethods(); document.body.classList.remove("sbopen"); }
    const typing = e.target.matches && e.target.matches("input,select,textarea");
    if ($("drawer").classList.contains("open") && !typing) {
      if (e.key === "ArrowDown" || e.key === "ArrowRight" || e.key === "j") { e.preventDefault(); stepDrawer(1); }
      else if (e.key === "ArrowUp" || e.key === "ArrowLeft" || e.key === "k") { e.preventDefault(); stepDrawer(-1); }
    }
    if (state.view !== "map" && !(e.target.matches && e.target.matches("input,select,textarea"))
        && !$("drawer").classList.contains("open")) {
      if (e.key === "ArrowRight") gotoPage(page + 1);
      else if (e.key === "ArrowLeft") gotoPage(page - 1);
      else if (e.key === "ArrowDown") { e.preventDefault(); setKbFocus(_kbIdx + 1); }
      else if (e.key === "ArrowUp") { e.preventDefault(); setKbFocus(_kbIdx - 1); }
      else if (e.key === "Enter") { const r = document.querySelectorAll("#tbody tr")[_kbIdx]; if (r) r.click(); }
    }
  });
  let _rz = null;
  addEventListener("resize", () => {
    setTopbarOffset();
    clearTimeout(_rz);
    _rz = setTimeout(() => { if (state.view === "table") renderThead(); }, 120);
  });
  loadOrg(state.org).then(() => {
    if (state._sel) { const r = DATA.rows.find((x) => x.uniprot_accession === state._sel); if (r) openDrawer(r); }
  });
}
init();
