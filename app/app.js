// app.js — target-prioritization browser (vanilla, no deps).
// Loads per-organism JSON, computes a live weighted composite, filters/sorts,
// renders a windowed OpenTargets-style table, and shows a per-protein detail drawer.
"use strict";

const LS_KEY = "gradi-tpb-state";
const PAGE_SIZE = 100; // rows per page

// min-threshold range filters (0–1). "__c" is the live composite.
const RANGES = [
  { key: "__c", label: "Composite ≥" },
  { key: "comp_essentiality", label: "Essentiality ≥" },
  { key: "comp_ligandability", label: "Ligandability ≥" },
  { key: "comp_breadth", label: "Breadth ≥" },
];

const state = {
  org: "kp",
  view: "overview",       // a TABLE_VIEWS key, or "map"
  mapColorBy: "__c",
  mapThreshold: 0,        // map: only highlight points with colour-by value >= this
  weights: {},            // key -> {on, weight}
  orthoTransfer: true,    // include ortholog-transferred knowledge (essentiality track)
  tiers: {},              // axisKey -> "any" | "low" | "med" | "high"
  filters: { search: "", ranges: {}, cats: {}, bools: {}, families: [] },
  shortlist: [],          // starred accessions (basket)
  shortlistOnly: false,   // filter table to the shortlist
  mapSel: null,           // Set of accessions box-selected on the map (or null)
  visibleCols: null,      // Set of column keys
  sort: { key: "__c", dir: -1 },
};

const cache = {};         // org -> parsed json
let DATA = null;          // current org payload
let AVAIL = {};           // component key -> bool (has data)
let filtered = [];        // current filtered+sorted rows
let page = 0;   // current page index (0-based)
let weightsDirty = true;  // recompute the composite only when weights actually change
let _lastRow = null;      // row currently shown in the drawer (for live re-render)

// ---------- persistence ----------------------------------------------------
function save() {
  const s = {
    org: state.org, view: state.view, mapColorBy: state.mapColorBy,
    tiers: state.tiers, mapThreshold: state.mapThreshold,
    weights: state.weights, orthoTransfer: state.orthoTransfer,
    filters: { search: state.filters.search, ranges: state.filters.ranges,
      cats: Object.fromEntries(Object.entries(state.filters.cats).map(([k, v]) => [k, [...v]])),
      bools: state.filters.bools, families: state.filters.families },
    shortlist: state.shortlist, shortlistOnly: state.shortlistOnly,
    visibleCols: state.visibleCols ? [...state.visibleCols] : null,
    sort: state.sort,
  };
  try { localStorage.setItem(LS_KEY, JSON.stringify(s)); } catch (e) {}
  syncHash();
}
function load() {
  try {
    const s = JSON.parse(localStorage.getItem(LS_KEY) || "null");
    if (!s) return;
    state.org = s.org || "kp";
    state.view = s.view || "overview";
    if (state.view === "table") state.view = "overview";   // migrate old value
    state.mapColorBy = s.mapColorBy || "__c";
    state.tiers = s.tiers || {};
    state.mapThreshold = typeof s.mapThreshold === "number" ? s.mapThreshold : 0;
    state.weights = s.weights || {};
    if (typeof s.orthoTransfer === "boolean") state.orthoTransfer = s.orthoTransfer;
    if (s.filters) {
      state.filters.search = s.filters.search || "";
      state.filters.ranges = s.filters.ranges || {};
      state.filters.cats = Object.fromEntries(Object.entries(s.filters.cats || {}).map(([k, v]) => [k, new Set(v)]));
      state.filters.bools = s.filters.bools || {};
      state.filters.families = Array.isArray(s.filters.families) ? s.filters.families : [];
    }
    state.shortlist = Array.isArray(s.shortlist) ? s.shortlist : [];
    state.shortlistOnly = !!s.shortlistOnly;
    state.visibleCols = s.visibleCols ? new Set(s.visibleCols) : null;
    state.sort = s.sort || state.sort;
  } catch (e) {}
}

// ---------- shareable deep-link (URL hash) ---------------------------------
let _applyingHash = false;
function serializeState() {
  const arr = (o) => Object.fromEntries(Object.entries(o).filter(([, v]) => v && v.length));
  return {
    org: state.org, view: state.view, xfer: state.orthoTransfer, map: state.mapColorBy,
    w: state.weights,
    t: arr(Object.fromEntries(Object.entries(state.tiers).map(([k, v]) => [k, Array.isArray(v) ? v : (v && v !== "any" ? [v] : [])]))),
    q: state.filters.search || undefined,
    cats: arr(Object.fromEntries(Object.entries(state.filters.cats).map(([k, v]) => [k, [...v]]))),
    bools: Object.fromEntries(Object.entries(state.filters.bools).filter(([, m]) => m && m !== "any")),
    fam: (state.filters.families && state.filters.families.length) ? state.filters.families : undefined,
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
    if (p.org) state.org = p.org;
    if (p.view) state.view = p.view;
    if (typeof p.xfer === "boolean") state.orthoTransfer = p.xfer;
    if (p.map) state.mapColorBy = p.map;
    if (p.w) state.weights = p.w;
    if (p.t) state.tiers = p.t;
    state.filters.search = p.q || "";
    state.filters.cats = Object.fromEntries(Object.entries(p.cats || {}).map(([k, v]) => [k, new Set(v)]));
    state.filters.bools = p.bools || {};
    state.filters.families = p.fam || [];
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
  const pfx = type === "sel" ? "sel-" : type === "pop" ? "pop-" : "tier-";
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
function enabledComponents() {
  return COMPONENTS.filter((c) => AVAIL[c.key] && state.weights[c.key] && state.weights[c.key].on
    && state.weights[c.key].weight > 0);
}
function computeAll() {
  const comps = enabledComponents();
  for (const row of DATA.rows) {
    let sum = 0, wsum = 0;
    for (const c of comps) {
      const v = row[c.key];
      if (isNum(v)) { const w = state.weights[c.key].weight; sum += w * v; wsum += w; }
    }
    row.__c = wsum > 0 ? sum / wsum : null;
    row.__conf = evidenceConfidence(row, comps).frac;
  }
}
// ---------- evidence / confidence ------------------------------------------
// Per-axis backing strength: 1 = measured/experimental, 0.5 = inferred
// (transfer / prediction / predicted pocket), 0 = present but unsupported
// (or provisional mock), null = axis has no value for this protein.
function axisBacking(row, key) {
  if (key === "comp_essentiality") {
    if (!isNum(row.comp_essentiality)) return null;
    return (row.experimentally_essential === true || isNum(row.evidence_experimental)) ? 1 : 0.5;
  }
  if (key === "comp_ligandability") {
    if (!isNum(row.comp_ligandability)) return null;
    return row.has_hard_evidence === true ? 1 : 0.5;
  }
  if (key === "comp_degradability") {
    if (!isNum(row.comp_degradability)) return null;
    if (isProvisional("comp_degradability")) return 0;         // mock data
    if (row.clp_trapped === true) return 1;
    return row.comp_degradability > 0 ? 0.5 : 0;
  }
  if (key === "comp_novelty") return isNum(row.comp_novelty) ? 1 : null;  // bibliometric = measured
  return isNum(row[key]) ? 0.5 : null;
}
function evidenceConfidence(row, comps) {
  comps = comps || enabledComponents();
  const per = []; let sum = 0, n = 0;
  for (const c of comps) {
    const b = axisBacking(row, c.key);
    if (b === null) continue;
    per.push({ label: c.label, backing: b }); sum += b; n++;
  }
  return { frac: n ? sum / n : null, per, n, predictedOnly: n > 0 && per.every((p) => p.backing < 1) };
}
function confGlyphHTML(row) {
  const ci = evidenceConfidence(row);
  if (ci.frac === null) return `<span class="conf" data-tip="No evidence for the enabled axes"><i></i><i></i><i></i></span>`;
  const col = ci.frac >= 0.8 ? "var(--good)" : ci.frac >= 0.4 ? "var(--warn)" : "var(--faint)";
  const dot = (t) => `<i class="${ci.frac >= t ? "f" : ""}"></i>`;
  const tip = ci.per.map((p) => `${p.label}: ${p.backing >= 1 ? "measured" : "inferred/predicted"}`).join(" · ")
    + (ci.predictedOnly ? " — predicted-only" : "");
  return `<span class="conf" style="--cf:${col}" data-tip="${tip.replace(/"/g, "&quot;")}">${dot(0.01)}${dot(0.5)}${dot(1)}</span>`;
}

// ---------- orthology-transfer toggle (essentiality axis) -----------------
// Essentiality evidence keys that come from cross-species transfer (dimmed in the
// card when transfer is off).
const TRANSFER_KEYS = new Set(["evidence_transfer", "ec_transfer_essential"]);
// Precompute, once per row: the transfer-included score (_ess_full, the stored value)
// and the transfer-excluded score (_ess_direct = renorm(0.40·experimental + 0.40·predictor)
// over present tracks — the exact essentiality_score formula minus the 0.20 transfer term).
function stashEssentiality() {
  for (const row of DATA.rows) {
    if (row._ess_full !== undefined) continue;
    row._ess_full = isNum(row.comp_essentiality) ? row.comp_essentiality : null;
    let num = 0, den = 0;
    if (isNum(row.evidence_experimental)) { num += 0.40 * row.evidence_experimental; den += 0.40; }
    if (isNum(row.evidence_predictor))   { num += 0.40 * row.evidence_predictor;   den += 0.40; }
    row._ess_direct = den > 0 ? num / den : null;
  }
}
function applyTransferMode() {
  for (const row of DATA.rows)
    row.comp_essentiality = state.orthoTransfer ? row._ess_full : row._ess_direct;
  weightsDirty = true;
}

// ---------- filtering / sorting -------------------------------------------
function passFilters(row) {
  const f = state.filters;
  if (state.shortlistOnly && !state.shortlist.includes(row.uniprot_accession)) return false;
  if (state.mapSel && !state.mapSel.has(row.uniprot_accession)) return false;
  if (f.search) {
    const q = f.search.toLowerCase();
    const hay = [row.gene, row.uniprot_accession, row.name, row.functional_class, row.pdb_ids,
      row.interpro_family_names, row.panther_family_names].filter(Boolean).join(" ").toLowerCase();
    if (!hay.includes(q)) return false;
  }
  for (const [key, min] of Object.entries(f.ranges)) {
    if (!min || min <= 0) continue;
    const v = row[key];
    if (!isNum(v) || v < min - 1e-9) return false;
  }
  for (const [key, set] of Object.entries(f.cats)) {
    if (!set || set.size === 0) continue;
    if (!set.has(row[key])) return false;
  }
  for (const [key, mode] of Object.entries(f.bools)) {
    if (!mode || mode === "any") continue;
    const want = mode === "yes";
    if (Boolean(row[key]) !== want) return false;
  }
  // family filter (OR across selected family names; matches any InterPro/PANTHER family field)
  if (f.families && f.families.length) {
    const fams = rowFamilySet(row);
    if (!f.families.some((name) => fams.has(name))) return false;
  }
  // tier selector (always active); multi-select = OR within an axis, AND across axes
  for (const axis of TIER_AXES) {
    const sel = tierSel(axis.key);
    if (!sel.length) continue;
    const v = row[axis.key];
    const bands = tierBands(axis);
    const inAny = sel.some((id) => {
      const b = bands.find((x) => x.id === id);
      return b && bandMatches(b, v, row);
    });
    if (!inAny) return false;
  }
  return true;
}
function recompute() {
  if (weightsDirty) { computeAll(); weightsDirty = false; }
  filtered = DATA.rows.filter(passFilters);
  const { key, dir } = state.sort;
  filtered.sort((a, b) => {
    let x = a[key], y = b[key];
    const xn = x === null || x === undefined || x === "", yn = y === null || y === undefined || y === "";
    if (xn && yn) return 0;
    if (xn) return 1;      // nulls always last
    if (yn) return -1;
    if (isNum(x) && isNum(y)) return (x - y) * dir;
    return String(x).localeCompare(String(y)) * dir;
  });
  page = 0;
  $("nShown").textContent = filtered.length.toLocaleString();
  $("nTotal").textContent = DATA.rows.length.toLocaleString();
  renderActiveView();
}
function renderActiveView() {
  const isMap = state.view === "map";
  $("tablewrap").hidden = isMap;
  $("mapwrap").hidden = !isMap;
  $("mapctrl").hidden = !isMap;
  $("colmenu").style.display = isMap ? "none" : "";
  $("pager").hidden = isMap;
  if (isMap) { drawMap(); const pn = $("provNote"); if (pn) pn.hidden = true; }
  else { renderThead(); renderPage(); }
  updateMapSelChip();
  setTopbarOffset();
}
function tableViewCols(key) {
  const v = TABLE_VIEWS.find((x) => x.key === key);
  if (!v) return [];
  return v.cols.filter((c) => DATA.columns.includes(c));
}
function buildViewSwitch() {
  const host = $("viewSwitch");
  const tabs = [...TABLE_VIEWS.map((v) => ({ view: v.key, label: v.label })),
                { view: "map", label: "◵ Map" }];
  host.innerHTML = tabs.map((t) =>
    `<button role="tab" data-view="${t.view}" aria-selected="${state.view === t.view}">${t.label}</button>`).join("");
  host.querySelectorAll("button").forEach((b) =>
    b.onclick = () => applyView(b.dataset.view));
}
function applyView(key) {
  state.view = key;
  $("viewSwitch").querySelectorAll("button").forEach((x) =>
    x.setAttribute("aria-selected", x.dataset.view === key));
  if (key !== "map") {
    state.visibleCols = new Set(tableViewCols(key));
    buildColMenu();
  }
  buildFilterBar();
  save();
  renderActiveView();
}

// ---------- table rendering (windowed) ------------------------------------
function visibleColumns() {
  const byKey = Object.fromEntries(TABLE_COLUMNS.map((c) => [c.key, c]));
  const ok = (k) => byKey[k] && DATA.columns.includes(k) && state.visibleCols.has(k);
  const view = TABLE_VIEWS.find((v) => v.key === state.view);
  const ordered = [];
  const seen = new Set();
  if (view) for (const k of view.cols) if (ok(k)) { ordered.push(byKey[k]); seen.add(k); }
  // append any extra columns the user toggled on that aren't in the view preset
  for (const c of TABLE_COLUMNS) if (!seen.has(c.key) && ok(c.key)) ordered.push(c);
  // keep same-section columns contiguous so a group band can never split in two
  // (e.g. a ligandability column toggled on via the Columns menu joins the
  // existing Ligandability block instead of forming a second band at the end).
  const groupsOrder = [];
  const byGroup = new Map();
  for (const c of ordered) {
    const g = columnGroup(c);
    if (!byGroup.has(g)) { byGroup.set(g, []); groupsOrder.push(g); }
    byGroup.get(g).push(c);
  }
  return groupsOrder.flatMap((g) => byGroup.get(g));
}
function tipAttr(text) {
  return text ? ` data-tip="${String(text).replace(/"/g, "&quot;")}"` : "";
}
const ROTATE_TYPES = new Set(["score", "pct", "plddt", "num", "int", "bool", "binary01"]);
function renderThead() {
  const cols = visibleColumns();
  const arrow = (k) => state.sort.key === k ? `<span class="arrow">${state.sort.dir < 0 ? "▼" : "▲"}</span>` : "";
  const th = (key, label, { nosort, tip, rot, color, cls: extra } = {}) => {
    const sortAttr = nosort ? "" : ` data-sort="${key}"`;
    const cls = [nosort ? "nosort" : "", rot ? "rot" : "", extra || ""].join(" ").trim();
    const styleAttr = rot ? ` style="--hc:${color || "var(--muted)"}"` : "";
    const inner = rot ? `<span class="rl">${label}${arrow(key)}</span>` : `${label}${arrow(key)}`;
    return `<th${sortAttr} class="${cls}"${styleAttr}${tipAttr(tip)}>${inner}</th>`;
  };
  // all columns are fixed-width with vertical headers, except the Target column
  let h = th("", "#", { nosort: true, tip: LEADING_COL_DESC.rank, cls: "colrank" })
    + th("name", "Target", { tip: LEADING_COL_DESC.name, cls: "colgene" })
    + th("__c", "Comp. score", { tip: LEADING_COL_DESC.__c, rot: true, color: AXIS_COLORS.composite })
    + th("__conf", "Evidence", { tip: LEADING_COL_DESC.__conf, rot: true, color: "var(--muted)" });
  for (const c of cols)
    h += th(c.key, c.label, { tip: c.desc, rot: true, color: colColor(c.key), cls: isProvisional(c.key) ? "provhdr" : "" });

  // section band: merge consecutive columns sharing a group
  const bandSeq = [{ key: "__c" }, { key: "__conf" }, ...cols];
  const runs = [];
  for (const col of bandSeq) {
    const g = columnGroup(col);
    const last = runs[runs.length - 1];
    if (last && last.group === g) last.span++;
    else runs.push({ group: g, span: 1 });
  }
  let band = `<th class="gband gband-empty" colspan="2"></th>`;
  for (const r of runs) {
    band += r.group
      ? `<th class="gband" colspan="${r.span}" style="--gc:${groupColor(r.group)}"><span>${r.group}</span></th>`
      : `<th class="gband" colspan="${r.span}"></th>`;
  }

  const thead = document.querySelector("#table thead");
  thead.innerHTML = `<tr class="bandrow">${band}</tr><tr>${h}</tr>`;
  thead.querySelectorAll("th[data-sort]").forEach((el) => {
    el.onclick = () => {
      const k = el.dataset.sort;
      if (state.sort.key === k) state.sort.dir *= -1;
      else state.sort = { key: k, dir: (k === "name") ? 1 : -1 };
      save(); recompute();
    };
  });
}
function tierType(key) {
  return key === "selectivity" ? "sel" : key === "popularity_tier" ? "pop" : "tier";
}
function classBadgeHTML(v) {
  if (!v) return "–";
  const m = FC_BY_ID[v];
  const col = m ? m.color : "var(--faint)";
  const abbr = m ? m.abbr : String(v).slice(0, 2);
  const label = m ? m.label : String(v).replace(/_/g, " ");
  return `<span class="badge" title="${label}" style="background:color-mix(in srgb,${col} 20%,#fff);color:color-mix(in srgb,${col} 75%,#000);">${abbr}</span>`;
}
function cellHTML(row, c) {
  const v = row[c.key];
  if (c.type === "bool") return boolHTML(v);
  if (c.type === "binary01") return (v === null || v === undefined) ? "–"
    : (v >= 0.5 ? '<span class="yes">✓</span>' : '<span class="cross">✗</span>');
  if (c.type === "class") return classBadgeHTML(v);
  if (c.type === "tier") return badgeHTML(tierType(c.key), v, c.key);
  if (c.type === "score" && c.heat) return `<span class="heat" style="${heatStyle(v, colColor(c.key))}">${fmt("score", v)}</span>`;
  if (c.type === "text") {
    if (v === null || v === undefined || v === "") return "–";
    const full = String(v).replace(/;/g, "; ");
    return `<span class="celltext" title="${full.replace(/"/g, "&quot;")}">${full}</span>`;
  }
  const cls = ["score", "pct", "plddt", "num", "int"].includes(c.type) ? "numcell" : "";
  return `<span class="${cls}">${fmt(c.type, v)}</span>`;
}
function pageCount() { return Math.max(1, Math.ceil(filtered.length / PAGE_SIZE)); }
function renderPage() {
  const cols = visibleColumns();
  if (page >= pageCount()) page = pageCount() - 1;
  if (page < 0) page = 0;
  const start = page * PAGE_SIZE;
  const end = Math.min(start + PAGE_SIZE, filtered.length);
  const frag = document.createDocumentFragment();
  for (let i = start; i < end; i++) {
    const row = filtered[i];
    const tr = document.createElement("tr");
    tr.dataset.acc = row.uniprot_accession;
    if (row.uniprot_accession === state._sel) tr.className = "sel";
    const cc = isNum(row.__c) ? row.__c : null;
    const fam = row.interpro_family_names || row.panther_family_names || row.interpro_superfamily_names || "";
    const fam1 = fam ? String(fam).split(";")[0] : "";
    const gtitle = fam1 ? ` title="${String(fam).replace(/;/g, "; ").replace(/"/g, "&quot;")}"` : "";
    let h = `<td class="colrank rank">${i + 1}</td>`
      + `<td class="colgene"${gtitle}>${starHTML(row.uniprot_accession)}<div class="gene">${row.name || row.gene || row.uniprot_accession}</div>`
      + `<div class="acc">${row.uniprot_accession}</div></td>`
      + `<td class="rc"><span class="heat cscore" style="${heatStyle(cc, AXIS_COLORS.composite)}">`
      + `${isNum(cc) ? cc.toFixed(2) : "–"}</span></td>`
      + `<td class="rc">${confGlyphHTML(row)}</td>`;
    for (const col of cols) h += `<td class="rc${isProvisional(col.key) ? " provcell" : ""}">${cellHTML(row, col)}</td>`;
    tr.innerHTML = h;
    tr.onclick = (e) => {
      const s = e.target.closest(".star");
      if (s) { e.stopPropagation(); toggleStar(row.uniprot_accession); } else openDrawer(row);
    };
    frag.appendChild(tr);
  }
  const tb = $("tbody");
  tb.innerHTML = "";
  tb.appendChild(frag);
  $("table").hidden = false;
  $("loading").hidden = true;
  window.scrollTo({ top: 0 });
  const pn = $("provNote");
  if (pn) pn.hidden = !cols.some((c) => isProvisional(c.key));
  renderPager();
  setTopbarOffset();
}
// ---------- methods / provenance modal -------------------------------------
function openMethods() {
  const m = $("methodsModal");
  m.innerHTML = `<div class="modalhead"><h2>How targets are scored</h2><button class="close" id="methodsClose" aria-label="Close">×</button></div>`
    + `<div class="modalbody">`
    + METHODS.map((s) => `<section class="mrow" style="--pc:${s.color}"><h3>${s.title}</h3><p>${s.body}</p></section>`).join("")
    + `</div>`;
  m.hidden = false; $("methodsScrim").hidden = false;
  $("methodsClose").onclick = closeMethods;
}
function closeMethods() { $("methodsModal").hidden = true; $("methodsScrim").hidden = true; }
// ---------- shortlist / basket ---------------------------------------------
function isStarred(acc) { return state.shortlist.includes(acc); }
function starHTML(acc) {
  return `<button class="star${isStarred(acc) ? " on" : ""}" data-acc="${acc}" title="Add to shortlist" aria-label="Shortlist">★</button>`;
}
function toggleStar(acc) {
  const i = state.shortlist.indexOf(acc);
  if (i >= 0) state.shortlist.splice(i, 1); else state.shortlist.push(acc);
  save(); updateShortlistBtn();
  if (state.shortlistOnly) recompute();
  else document.querySelectorAll(`.star[data-acc="${acc}"]`).forEach((el) => el.classList.toggle("on", isStarred(acc)));
}
function updateShortlistBtn() {
  const b = $("shortlistBtn");
  if (!b) return;
  b.textContent = `★ Shortlist (${state.shortlist.length})`;
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
}
function renderPager() {
  const host = $("pager");
  const n = pageCount();
  host.hidden = state.view === "map";
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
// SVG donut for a 0–1 score (arc only; the number is overlaid in HTML)
function ringSVG(v, color) {
  const r = 42, circ = 2 * Math.PI * r;
  const val = isNum(v) ? Math.max(0, Math.min(1, v)) : 0;
  const off = circ * (1 - val);
  return `<svg viewBox="0 0 100 100" class="ringsvg" aria-hidden="true">`
    + `<circle cx="50" cy="50" r="${r}" class="ringbg"></circle>`
    + `<circle cx="50" cy="50" r="${r}" class="ringfg" style="stroke:${color};`
    + `stroke-dasharray:${circ.toFixed(1)};stroke-dashoffset:${off.toFixed(1)}"></circle></svg>`;
}
// AlphaFold pLDDT confidence palette
function plddtColor(v) {
  if (!isNum(v)) return "var(--faint)";
  if (v >= 90) return "#0053D6"; if (v >= 70) return "#65CBF3";
  if (v >= 50) return "#FFDB13"; return "#FF7D45";
}
function plddtBarHTML(v) {
  const pct = Math.max(0, Math.min(100, v));
  return `<div class="plddt"><div class="plddtscale"></div>`
    + `<div class="plddtneedle" style="left:${pct}%" title="mean pLDDT ${Math.round(v)}"></div>`
    + `<div class="plddtlabels"><span>0</span><span>disorder → confident</span><span>100</span></div></div>`;
}
// labelled 0–1 metric bar (invert = higher is worse → risk colour)
function barHTML(label, v, color, invert, small) {
  const has = isNum(v);
  const pct = has ? Math.round(Math.max(0, Math.min(1, v)) * 100) : 0;
  const bcol = invert ? "var(--tangerine)" : color;
  return `<div class="metric${small ? " mini" : ""}"><span class="ml">${label}</span>`
    + `<span class="bar"><i style="width:${pct}%;background:${bcol}"></i></span>`
    + `<span class="mv">${has ? v.toFixed(2) : "–"}</span></div>`;
}
function statHTML(label, v, type) {
  const shown = (v === null || v === undefined || v === "") ? "–" : fmt(type, v);
  return `<div class="stat2"><div class="s2v">${esc(shown)}</div><div class="s2l">${label}</div></div>`;
}
function flagHTML(label, on) {
  return `<span class="flag ${on ? "on" : "off"}">${on ? "✓" : "·"} ${label}</span>`;
}
function pdbChipsHTML(str) {
  const ids = String(str || "").split(/[;,\s]+/).filter(Boolean).slice(0, 12);
  return ids.map((id) => `<a class="chip pdbchip" href="https://www.rcsb.org/structure/${encodeURIComponent(id)}" target="_blank" rel="noopener">${esc(id)}</a>`).join("");
}
function chipClusterHTML(str, max) {
  const parts = String(str || "").split(";").map((s) => s.trim()).filter(Boolean);
  if (!parts.length) return "";
  const shown = parts.slice(0, max || 10);
  const extra = parts.length > shown.length ? `<span class="chip more">+${parts.length - shown.length}</span>` : "";
  return shown.map((p) => `<span class="chip">${esc(p)}</span>`).join("") + extra;
}
function tierChipHTML(type, v) {
  if (v === null || v === undefined || v === "") return "";
  const pfx = type === "sel" ? "sel-" : type === "pop" ? "pop-" : "tier-";
  return `<span class="badge ${pfx}${v}">${esc(String(v).replace(/_/g, " "))}</span>`;
}
function axisPanelHTML(row, spec) {
  const has = (k) => DATA.columns.includes(k);
  const color = CARD_GROUP_COLORS[spec.axis] || "var(--brand)";
  // header: title + tier badge + headline score
  let meta = "";
  if (spec.tier && has(spec.tier) && row[spec.tier] != null) meta += tierChipHTML(tierType(spec.tier), row[spec.tier]);
  else if (spec.tierClass && has(spec.tierClass) && row[spec.tierClass] != null) meta += tierChipHTML("sel", row[spec.tierClass]);
  if (spec.headline && has(spec.headline)) {
    const v = row[spec.headline];
    meta += `<span class="pscore" style="color:${color}">${isNum(v) ? v.toFixed(2) : "–"}</span>`;
  } else if (spec.plddt && has(spec.plddt)) {
    const v = row[spec.plddt];
    meta += `<span class="pscore" style="color:${plddtColor(v)}">${isNum(v) ? Math.round(v) : "–"}<small>pLDDT</small></span>`;
  }
  const prov = (spec.provisionalOrgs && spec.provisionalOrgs.includes(state.org))
    ? `<span class="provchip">provisional · mock</span>` : "";
  const head = `<div class="phead"><h3 style="--pc:${color}">${spec.title}</h3><div class="pmeta">${prov}${meta}</div></div>`;
  // dim cross-species transfer evidence when the orthology-transfer toggle is off
  const excluded = (k) => !state.orthoTransfer && TRANSFER_KEYS.has(k);

  let body = "";
  if (spec.blurb) body += `<p class="pblurb">${spec.blurb}</p>`;
  if (spec.bars) for (const [k, label, invert] of spec.bars) if (has(k)) {
    let h = barHTML(label, row[k], color, invert);
    if (excluded(k)) h = h.replace('class="metric"', 'class="metric excluded"');
    body += h;
  }
  if (spec.subbars) {
    const sb = spec.subbars.filter(([k]) => has(k) && row[k] != null);
    if (sb.length) body += `<div class="subbars">` + sb.map(([k, label]) => barHTML(label, row[k], color, false, true)).join("") + `</div>`;
  }
  if (spec.plddt && has(spec.plddt) && isNum(row[spec.plddt])) body += plddtBarHTML(row[spec.plddt]);
  if (spec.penalty && has(spec.penalty[0]) && isNum(row[spec.penalty[0]])) body += barHTML(spec.penalty[1], row[spec.penalty[0]], "var(--silver)", true);
  if (spec.stats) {
    const st = spec.stats.filter(([k]) => has(k));
    if (st.length) body += `<div class="statgrid">` + st.map(([k, label, type]) => statHTML(label, row[k], type)).join("") + `</div>`;
  }
  if (spec.flags) {
    const fl = spec.flags.filter(([k]) => has(k));
    if (fl.length) body += `<div class="flagrow">` + fl.map(([k, label, mode]) => {
      const raw = row[k];
      const on = mode === "ge05" ? (isNum(raw) && raw >= 0.5) : Boolean(raw);
      const f = flagHTML(label, on);
      return excluded(k) ? f.replace('class="flag ', 'class="flag excluded ') : f;
    }).join("") + `</div>`;
  }
  if (spec.chips) for (const [k, label, kind] of spec.chips) {
    if (!has(k) || row[k] === null || row[k] === undefined || row[k] === "") continue;
    const inner = kind === "pdb" ? pdbChipsHTML(row[k]) : `<span class="chip">${esc(String(row[k]))}</span>`;
    if (inner) body += `<div class="chiprow"><span class="crl">${label}</span><span class="chips">${inner}</span></div>`;
  }
  if (spec.text) for (const [k, label] of spec.text) {
    if (!has(k) || row[k] === null || row[k] === undefined || row[k] === "") continue;
    body += `<div class="kv"><span class="kvk">${label}</span><span class="kvv">${esc(String(row[k]))}</span></div>`;
  }
  if (spec.homolog) {
    const g = row[spec.homolog[0]], o = row[spec.homolog[1]];
    if ((g && g !== "") || (o && o !== "")) body += `<div class="kv"><span class="kvk">Best-studied homolog</span>`
      + `<span class="kvv">${esc(g || "–")}${o ? ` <em>${esc(o)}</em>` : ""}</span></div>`;
  }
  if (spec.sources && has(spec.sources) && row[spec.sources]) {
    const chips = chipClusterHTML(row[spec.sources], 12);
    if (chips) body += `<div class="chiprow"><span class="crl">Evidence</span><span class="chips">${chips}</span></div>`;
  }
  if (!body.trim()) return "";
  return `<section class="panel">${head}${body}</section>`;
}
function annotationPanelHTML(row) {
  const has = (k) => DATA.columns.includes(k);
  const spec = CARD_ANNOTATION, color = CARD_GROUP_COLORS.annotation;
  let body = "";
  const fc = row[spec.classKey];
  if (fc) {
    const m = FC_BY_ID[fc];
    body += `<div class="kv"><span class="kvk">Functional class</span>`
      + `<span class="kvv">${classBadgeHTML(fc)} ${m ? m.label : esc(String(fc))}</span></div>`;
  }
  for (const [k, label] of spec.chipGroups) {
    if (!has(k) || !row[k]) continue;
    const chips = chipClusterHTML(row[k]);
    if (chips) body += `<div class="chiprow"><span class="crl">${label}</span><span class="chips">${chips}</span></div>`;
  }
  if (!body.trim()) return "";
  return `<section class="panel"><div class="phead"><h3 style="--pc:${color}">${spec.title}</h3></div>${body}</section>`;
}
// mini ESM-C locator map with this protein highlighted
function drawCardMap(row) {
  const cv = document.getElementById("cardmap");
  if (!cv || !DATA) return;
  const W = cv.clientWidth || 460, H = cv.clientHeight || 172, dpr = devicePixelRatio || 1;
  cv.width = W * dpr; cv.height = H * dpr;
  const ctx = cv.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);
  const { minx, maxx, miny, maxy } = mapExtent();
  const P = 12, s = Math.min((W - 2 * P) / ((maxx - minx) || 1), (H - 2 * P) / ((maxy - miny) || 1));
  const ox = P + (W - 2 * P - s * (maxx - minx)) / 2, oy = P + (H - 2 * P - s * (maxy - miny)) / 2;
  const X = (x) => ox + (x - minx) * s, Y = (y) => oy + (maxy - y) * s;
  ctx.fillStyle = "rgba(154,147,166,0.20)";
  for (const r of DATA.rows) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    ctx.beginPath(); ctx.arc(X(r.tsne_x), Y(r.tsne_y), 1, 0, 6.283); ctx.fill();
  }
  if (isNum(row.tsne_x) && isNum(row.tsne_y)) {
    const px = X(row.tsne_x), py = Y(row.tsne_y);
    ctx.strokeStyle = "rgba(108,92,231,0.30)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(px, P); ctx.lineTo(px, H - P); ctx.moveTo(P, py); ctx.lineTo(W - P, py); ctx.stroke();
    ctx.beginPath(); ctx.arc(px, py, 5.5, 0, 6.283);
    ctx.fillStyle = "#6C5CE7"; ctx.fill();
    ctx.lineWidth = 2; ctx.strokeStyle = "#fff"; ctx.stroke();
  }
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
async function loadPockets(org) {
  if (pocketCache[org]) return pocketCache[org];
  try {
    const r = await fetch(`data/pockets_${org}.json`);
    pocketCache[org] = r.ok ? await r.json() : {};
  } catch (e) { pocketCache[org] = {}; }
  return pocketCache[org];
}
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
  if (row.af_available === false) {
    msg(`No AlphaFold model for this protein.`);
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
      v.zoomTo(sel); v.zoom(0.85);
    } else {
      v.zoomTo();
    }
    v.render();
    const cap = document.getElementById("pdbcap");
    if (cap) cap.innerHTML = (resids && resids.length)
      ? `AlphaFold model · cartoon coloured by pLDDT · <span class="pk">green = top predicted pocket</span> (${resids.length} residues). Drag to rotate, scroll to zoom.`
      : `AlphaFold model · cartoon coloured by pLDDT. Drag to rotate, scroll to zoom.`;
  };
  requestAnimationFrame(paint);
}
function openDrawer(row) {
  state._sel = row.uniprot_accession;
  _lastRow = row;
  document.querySelectorAll("#tbody tr").forEach((tr) =>
    tr.classList.toggle("sel", tr.dataset.acc === row.uniprot_accession));

  const has = (k) => DATA.columns.includes(k);
  const cval = isNum(row.__c) ? row.__c : 0;
  const rank = compRank(row);

  // composite contribution bars (per enabled component)
  const comps = COMPONENTS.filter((c) => AVAIL[c.key] && state.weights[c.key] && state.weights[c.key].on
    && state.weights[c.key].weight > 0);
  const contribHTML = comps.map((c) => {
    const v = row[c.key];
    return `<div class="c" style="--cc:${colColor(c.key)}"><span class="cl">${c.label}</span>`
      + `<div class="cbarmini"><i style="width:${isNum(v) ? Math.round(v * 100) : 0}%"></i></div>`
      + `<span class="cvv">${isNum(v) ? v.toFixed(2) : "–"}</span></div>`;
  }).join("") || '<span class="dk">No active components</span>';

  // functional-class pill
  const fcMeta = FC_BY_ID[row.functional_class];
  const fcPill = fcMeta ? `<span class="fcpill" style="--fc:${fcMeta.color}" title="Functional class">${fcMeta.label}</span>` : "";

  // at-a-glance tier chips
  const glance = [
    tierChipHTML("tier", has("essentiality_tier") ? row.essentiality_tier : null),
    tierChipHTML("tier", has("ligandability_tier") ? row.ligandability_tier : null),
    tierChipHTML("sel", has("selectivity") ? row.selectivity : null),
    tierChipHTML("pop", has("popularity_tier") ? row.popularity_tier : null),
  ].filter(Boolean).join("");
  const glanceHTML = glance ? `<div class="tierchips">${glance}</div>` : "";

  const links = externalLinks(row).map((l) =>
    `<a href="${l.href}" target="_blank" rel="noopener">${l.icon ? `<span class="li">${l.icon}</span>` : ""}${l.label}</a>`).join("");

  const panels = CARD_AXES.map((spec) => axisPanelHTML(row, spec)).join("");
  const showViewer = !has("af_available") || row.af_available !== false;
  const viewerPanel = showViewer
    ? `<section class="panel viewerpanel"><div class="phead"><h3 style="--pc:${CARD_GROUP_COLORS.ligandability}">3D structure</h3></div>`
      + `<div class="pdbslot" id="pdbslot"><div class="pdbmsg">Loading 3D structure…</div></div>`
      + `<p class="pdbcap" id="pdbcap"></p></section>`
    : "";
  const mapPanel = (has("tsne_x") && isNum(row.tsne_x))
    ? `<section class="panel minipanel"><div class="phead"><h3 style="--pc:var(--brand)">Protein universe</h3></div>`
      + `<p class="pblurb">Where this protein sits in the ESM-C embedding map (all ${DATA.rows.length.toLocaleString()} proteins).</p>`
      + `<div class="minimap"><canvas id="cardmap"></canvas></div></section>`
    : "";
  const annotation = annotationPanelHTML(row);

  $("drawer").innerHTML = `
    <div class="dhead">
      <button class="close" id="drawerClose" aria-label="Close">×</button>
      <div class="idrow">${starHTML(row.uniprot_accession)}<span class="gene">${esc(row.gene || row.uniprot_accession)}</span>${fcPill}</div>
      <div class="submeta">
        <button class="accbtn" id="accCopy" title="Copy accession">${esc(row.uniprot_accession)} <span class="cpi">⧉</span></button>
        <span class="org"><span class="odot odot-${state.org}"></span>${ORGANISM_META[state.org].name} ${ORGANISM_META[state.org].strain}</span>
      </div>
      <div class="scorehero">
        <div class="ring">${ringSVG(cval, "var(--brand)")}<div class="ringc"><b>${isNum(row.__c) ? row.__c.toFixed(2) : "–"}</b><span>composite</span></div></div>
        <div class="herostack">
          ${rank ? `<div class="rankline">Rank <b>#${rank.toLocaleString()}</b> <span>of ${DATA.rows.length.toLocaleString()}</span>`
            + `<span class="evline">${confGlyphHTML(row)} ${(function(){const ci=evidenceConfidence(row);return ci.predictedOnly?'<b class="predonly">predicted-only</b>':(ci.frac>=0.8?'well supported':ci.frac>=0.4?'partly supported':'weakly supported');})()}</span></div>` : ""}
          <div class="contrib">${contribHTML}</div>
        </div>
      </div>
      ${glanceHTML}
      <div class="links">${links}</div>
    </div>
    <div class="dbody">${viewerPanel}${panels}${mapPanel}${annotation}</div>`;
  $("drawerClose").onclick = closeDrawer;
  const cardStar = $("drawer").querySelector(".idrow .star");
  if (cardStar) cardStar.onclick = () => toggleStar(row.uniprot_accession);
  const accBtn = $("accCopy");
  if (accBtn) accBtn.onclick = () => {
    try { navigator.clipboard && navigator.clipboard.writeText(row.uniprot_accession); } catch (e) {}
    accBtn.classList.add("copied"); setTimeout(() => accBtn.classList.remove("copied"), 1100);
  };
  $("drawer").classList.add("open");
  $("drawer").setAttribute("aria-hidden", "false");
  $("scrim").classList.add("open");
  requestAnimationFrame(() => drawCardMap(row));
  if (showViewer) renderStructure(row);
}
function closeDrawer() {
  $("drawer").classList.remove("open");
  $("drawer").setAttribute("aria-hidden", "true");
  $("scrim").classList.remove("open");
  state._sel = null;
  document.querySelectorAll("#tbody tr.sel").forEach((tr) => tr.classList.remove("sel"));
}

// ---------- selection mode (weights vs tiers) -----------------------------
// normalize a tier selection to an array of band ids (migrates old single-string form)
function tierSel(axisKey) {
  const v = state.tiers[axisKey];
  if (Array.isArray(v)) return v;
  if (typeof v === "string" && v !== "any") return [v];
  return [];
}
function tierCount(axisKey, bandId) {
  const axis = TIER_AXES.find((a) => a.key === axisKey);
  const band = tierBands(axis).find((b) => b.id === bandId);
  let n = 0;
  for (const r of DATA.rows) if (bandMatches(band, r[axisKey], r)) n++;
  return n;
}
function buildTierPanel() {
  const host = $("tierPanel");
  host.innerHTML = "";
  for (const axis of TIER_AXES) {
    if (!DATA.columns.includes(axis.key)) continue;
    const sel = tierSel(axis.key);          // array of selected band ids (multi-select)
    const col = colColor(axis.key);
    const div = document.createElement("div");
    div.className = "tieraxis";
    div.style.setProperty("--tc", col);
    const bands = tierBands(axis);
    const hasSel = sel.length > 0;
    const segs = bands.map((b) =>
      `<button data-band="${b.id}" title="${b.label} (${tierCount(axis.key, b.id).toLocaleString()} proteins)" `
      + `aria-selected="${sel.includes(b.id)}" style="--bi:${b.intensity}%">${b.label}</button>`).join("");
    div.innerHTML = `
      <div class="tl"><span class="dot"></span><span class="tlabel">${axis.label}</span>
        <button class="tany" data-band="any" aria-selected="${!hasSel}">any</button></div>
      <div class="bands${hasSel ? " hassel" : ""}" style="--n:${bands.length}">${segs}</div>`;
    div.querySelectorAll("button[data-band]").forEach((b) => {
      b.onclick = () => {
        const cur = tierSel(axis.key);
        if (b.dataset.band === "any") {
          state.tiers[axis.key] = [];
        } else if (cur.includes(b.dataset.band)) {
          state.tiers[axis.key] = cur.filter((x) => x !== b.dataset.band);
        } else {
          state.tiers[axis.key] = [...cur, b.dataset.band];
        }
        buildTierPanel(); save(); recompute();
      };
    });
    host.appendChild(div);
  }
}

// ---------- sidebar builders ----------------------------------------------
function buildWeights() {
  const host = $("weights");
  host.innerHTML = "";
  for (const c of COMPONENTS) {
    const avail = AVAIL[c.key];
    if (!state.weights[c.key]) state.weights[c.key] = { on: c.on, weight: c.weight };
    const w = state.weights[c.key];
    const div = document.createElement("div");
    div.className = "weight" + (!avail ? " disabled" : (w.on ? "" : " off"));
    div.style.setProperty("--wc", colColor(c.key));
    div.innerHTML = `
      <div class="top">
        <input type="checkbox" class="sw" ${w.on && avail ? "checked" : ""} ${avail ? "" : "disabled"} aria-label="Enable ${c.label} in composite">
        <span class="lab">${c.label}${c.binary ? ` <span class="btag" title="Binary axis: value is 0 or 1. The weight is added only to targets scoring 1.">0/1</span>` : ""}</span>
        ${avail ? `<span class="wv">${w.weight}</span>` : `<span class="soon">soon</span>`}
      </div>
      ${avail ? `<input type="range" min="0" max="100" step="5" value="${w.weight}" aria-label="${c.label} weight">` : ""}
      <div class="help">${c.help}</div>`;
    if (avail) {
      const sw = div.querySelector(".sw"), rng = div.querySelector("input[type=range]"), wv = div.querySelector(".wv");
      sw.onchange = () => { w.on = sw.checked; div.classList.toggle("off", !w.on); weightsDirty = true; save(); recompute(); };
      rng.oninput = () => { w.weight = +rng.value; wv.textContent = rng.value; };
      rng.onchange = () => { weightsDirty = true; save(); recompute(); };
    }
    host.appendChild(div);
  }
}
function distinct(key) {
  const s = new Set();
  for (const r of DATA.rows) if (r[key] !== null && r[key] !== undefined && r[key] !== "") s.add(r[key]);
  return [...s].sort();
}
// One categorical chip group as an inline filter-bar element.
function catFilterEl(key) {
  const f = CATEGORICAL_FILTERS.find((x) => x.key === key) || { key, label: key };
  if (!state.filters.cats[key]) state.filters.cats[key] = new Set();
  const active = state.filters.cats[key];
  const div = document.createElement("div");
  div.className = "fgroup";
  div.innerHTML = `<span class="flabel">${f.label}</span><div class="chk"></div>`;
  const chk = div.querySelector(".chk");
  const present = new Set(distinct(key));
  const values = f.order ? f.order.filter((v) => present.has(v)) : [...present];
  for (const val of values) {
    const meta = f.meta && f.meta[val];
    const b = document.createElement("button");
    b.textContent = meta ? meta.label : String(val).replace(/_/g, " ");
    b.setAttribute("aria-pressed", active.has(val) ? "true" : "false");
    if (meta) { b.style.setProperty("--cc", meta.color); b.classList.add("cchip"); }
    b.onclick = () => {
      if (active.has(val)) active.delete(val); else active.add(val);
      b.setAttribute("aria-pressed", active.has(val) ? "true" : "false");
      save(); recompute();
    };
    chk.appendChild(b);
  }
  return div;
}
// One tri-state boolean as an inline filter-bar element.
function boolFilterEl(key) {
  const f = BOOL_FILTERS.find((x) => x.key === key) || { key, label: key };
  if (!state.filters.bools[key]) state.filters.bools[key] = "any";
  const mode = state.filters.bools[key];
  const div = document.createElement("div");
  div.className = "fgroup bool-filter";
  div.innerHTML = `<span class="flabel">${f.label}</span><div class="tristate">`
    + ["any", "yes", "no"].map((m) => `<button data-m="${m}" aria-pressed="${mode === m}">${m}</button>`).join("")
    + `</div>`;
  div.querySelectorAll("button").forEach((b) => {
    b.onclick = () => {
      state.filters.bools[key] = b.dataset.m;
      div.querySelectorAll("button").forEach((x) => x.setAttribute("aria-pressed", x === b));
      save(); recompute();
    };
  });
  return div;
}
// family names live across several ;-joined columns
const FAMILY_FIELDS = ["interpro_family_names", "interpro_superfamily_names", "panther_family_names"];
function rowFamilySet(row) {
  const s = new Set();
  for (const k of FAMILY_FIELDS) {
    const v = row[k];
    if (typeof v === "string" && v) for (const p of v.split(";")) { const t = p.trim(); if (t) s.add(t); }
  }
  return s;
}
let FAMILY_VOCAB = [];   // distinct family names, sorted by frequency (for autocomplete)
function buildFamilyVocab() {
  const counts = new Map();
  for (const r of DATA.rows) for (const name of rowFamilySet(r)) counts.set(name, (counts.get(name) || 0) + 1);
  FAMILY_VOCAB = [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([n]) => n);
}
// Protein-family autocomplete as an inline filter-bar element.
function familyFilterEl() {
  const sel = state.filters.families;
  const opts = FAMILY_VOCAB.slice(0, 1200).map((n) => `<option value="${n.replace(/"/g, "&quot;")}">`).join("");
  const div = document.createElement("div");
  div.className = "fgroup famgroup";
  div.innerHTML = `<span class="flabel">Protein family</span>`
    + `<input type="text" class="famInput" list="familyList" placeholder="type to search ${FAMILY_VOCAB.length.toLocaleString()}…" autocomplete="off" aria-label="Protein family">`
    + `<datalist id="familyList">${opts}</datalist>`
    + `<div class="chk famchips"></div>`;
  const chips = div.querySelector(".famchips");
  const renderChips = () => {
    chips.innerHTML = sel.map((n, i) =>
      `<button data-i="${i}" title="${n.replace(/"/g, "&quot;")}" aria-pressed="true">${n} ✕</button>`).join("");
    chips.querySelectorAll("button").forEach((b) => {
      b.onclick = () => { sel.splice(+b.dataset.i, 1); renderChips(); save(); recompute(); };
    });
  };
  const inp = div.querySelector(".famInput");
  const add = () => {
    const val = inp.value.trim();
    if (val && FAMILY_VOCAB.includes(val) && !sel.includes(val)) {
      sel.push(val); inp.value = ""; renderChips(); save(); recompute();
    }
  };
  inp.onchange = add;
  inp.onkeydown = (e) => { if (e.key === "Enter") { e.preventDefault(); add(); } };
  renderChips();
  return div;
}
function presetsEl() {
  const div = document.createElement("div");
  div.className = "fgroup presets";
  div.innerHTML = `<button data-preset="prime">★ Prime targets</button>`
    + `<button data-preset="neglected">◐ Neglected &amp; druggable</button>`;
  div.querySelectorAll("button").forEach((b) => b.onclick = () => applyPreset(b.dataset.preset));
  return div;
}
function anyFilterActive() {
  const f = state.filters;
  return !!(f.search || (f.families && f.families.length)
    || Object.values(f.cats).some((s) => s && s.size)
    || Object.values(f.bools).some((m) => m && m !== "any")
    || TIER_AXES.some((a) => tierSel(a.key).length));
}
function clearEl() {
  const div = document.createElement("div");
  div.className = "fgroup";
  const b = document.createElement("button");
  b.className = "clearbtn"; b.textContent = "Clear filters";
  b.onclick = () => applyPreset("reset");
  div.appendChild(b);
  return div;
}
// Render the per-page filter bar (main area) for the active view.
function buildFilterBar() {
  const host = $("filterbar");
  const spec = VIEW_FILTERS[state.view] || {};
  host.innerHTML = "";
  if (state.view === "map") { host.hidden = true; return; }
  const parts = [];
  if (spec.presets) parts.push(presetsEl());
  (spec.cats || []).forEach((k) => { if (DATA.columns.includes(k)) parts.push(catFilterEl(k)); });
  if (spec.family && DATA.columns.some((c) => FAMILY_FIELDS.includes(c))) parts.push(familyFilterEl());
  (spec.bools || []).forEach((k) => { if (DATA.columns.includes(k)) parts.push(boolFilterEl(k)); });
  if (parts.length) parts.push(clearEl());
  parts.forEach((p) => host.appendChild(p));
  host.hidden = parts.length === 0;
}
function buildColMenu() {
  const pop = $("colPop");
  pop.innerHTML = "";
  for (const c of TABLE_COLUMNS) {
    if (!DATA.columns.includes(c.key)) continue;
    const lab = document.createElement("label");
    const on = state.visibleCols.has(c.key);
    lab.innerHTML = `<input type="checkbox" ${on ? "checked" : ""}> ${c.label}`;
    lab.querySelector("input").onchange = (e) => {
      if (e.target.checked) state.visibleCols.add(c.key); else state.visibleCols.delete(c.key);
      save(); renderThead(); renderPage();
    };
    pop.appendChild(lab);
  }
}

// ---------- presets --------------------------------------------------------
function applyPreset(name) {
  if (name === "reset") {
    state.filters.search = ""; if ($("search")) $("search").value = "";
    state.filters.ranges = {};
    for (const k in state.filters.cats) state.filters.cats[k].clear();
    for (const k in state.filters.bools) state.filters.bools[k] = "any";
    state.filters.families = [];
    for (const a of TIER_AXES) state.tiers[a.key] = [];
    buildTierPanel();
  } else if (name === "prime") {
    applyPreset("reset");
    state.filters.cats.essentiality_tier = new Set(["essential"]);
    state.filters.cats.ligandability_tier = new Set(["tractable"]);
    state.filters.cats.selectivity = new Set(["broad_selective"]);
  } else if (name === "neglected") {
    applyPreset("reset");
    // druggable + understudied: tractable ligandability crossed with dark studiedness
    state.filters.cats.ligandability_tier = new Set(["tractable"]);
    if (DATA.columns.includes("popularity_tier")) state.filters.cats.popularity_tier = new Set(["dark"]);
  }
  buildFilterBar();
  save(); recompute();
}

// ---------- CSV export -----------------------------------------------------
function exportCSV() {
  const cols = visibleColumns();
  const header = ["rank", "uniprot_accession", "gene", "composite", "evidence_support", ...cols.map((c) => c.key)];
  const lines = [header.join(",")];
  filtered.forEach((row, i) => {
    const vals = [i + 1, row.uniprot_accession, row.gene || "", isNum(row.__c) ? row.__c.toFixed(4) : "",
      isNum(row.__conf) ? row.__conf.toFixed(2) : ""];
    for (const c of cols) { let v = row[c.key]; v = (v === null || v === undefined) ? "" : v; vals.push(v); }
    lines.push(vals.map((v) => {
      const s = String(v);
      return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
    }).join(","));
  });
  const blob = new Blob([lines.join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `gradi_${state.org}_targets.csv`;
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
async function loadOrg(org) {
  state.org = org;
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.setAttribute("aria-selected", b.dataset.org === org));
  $("table").hidden = true; $("loading").hidden = false; $("loading").textContent = "Loading…";
  if (!cache[org]) {
    const meta = ORGANISM_META[org];
    const res = await fetch(meta.file);
    if (!res.ok) { $("loading").textContent = `Failed to load ${meta.file} (${res.status})`; return; }
    cache[org] = rehydrate(await res.json());
  }
  DATA = cache[org];
  weightsDirty = true;   // composite must be (re)computed for the newly-active organism
  AVAIL = {};
  for (const c of COMPONENTS) AVAIL[c.key] = DATA.columns.includes(c.key) && (c.available !== false);
  // seed visible columns from the active table view (fall back to overview)
  const seedView = state.view === "map" ? "overview" : state.view;
  if (!state.visibleCols || state.visibleCols.size === 0) state.visibleCols = new Set(tableViewCols(seedView));
  $("genInfo").innerHTML = `Data generated <code>${DATA.generated_at}</code> · ${DATA.n.toLocaleString()} proteins · `
    + `${ORGANISM_META[org].name} ${ORGANISM_META[org].strain}`;
  buildViewSwitch();
  buildWeights(); buildTierPanel();
  buildFamilyVocab();
  buildFilterBar(); buildColMenu();
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

  const t = state.mapThreshold || 0;
  // base layer: every protein as faint context; also collect hit-test points
  mapPts = [];
  ctx.fillStyle = `rgba(${faint.r},${faint.g},${faint.b},0.16)`;
  for (const r of DATA.rows) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    const px = X(r.tsne_x), py = Y(r.tsne_y);
    ctx.beginPath(); ctx.arc(px, py, 1.5, 0, 6.283); ctx.fill();
    mapPts.push({ sx: px, sy: py, row: r });
  }
  // highlight layer: passing rows whose colour-by value clears the threshold
  let hi = 0;
  ctx.lineWidth = 0.6; ctx.strokeStyle = "rgba(255,255,255,0.85)";
  for (const r of filtered) {
    if (!isNum(r.tsne_x) || !isNum(r.tsne_y)) continue;
    const v = r[colorBy];
    if (!isNum(v) || v < t) continue;
    const vv = Math.max(0, Math.min(1, v));
    ctx.fillStyle = `rgba(${acc.r},${acc.g},${acc.b},${0.5 + 0.5 * vv})`;
    ctx.beginPath(); ctx.arc(X(r.tsne_x), Y(r.tsne_y), 3 + 1.8 * vv, 0, 6.283); ctx.fill(); ctx.stroke();
    hi++;
  }
  const cnt = $("mapHiCount"); if (cnt) cnt.textContent = `${hi.toLocaleString()} highlighted`;
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
  if (state._sel) mark(state._sel, 7);
  if (_hoverAcc && _hoverAcc !== state._sel) mark(_hoverAcc, 6);
  if (_mapDrag) {
    const { x0, y0, x1, y1 } = _mapDrag;
    ctx.fillStyle = "rgba(108,92,231,0.10)"; ctx.strokeStyle = "rgba(108,92,231,0.9)"; ctx.lineWidth = 1;
    const rx = Math.min(x0, x1), ry = Math.min(y0, y1), rw = Math.abs(x1 - x0), rh = Math.abs(y1 - y0);
    ctx.fillRect(rx, ry, rw, rh); ctx.strokeRect(rx, ry, rw, rh);
  }
}
function renderMapLegend() {
  const lbl = (MAP_COLORS.find((m) => m.key === state.mapColorBy) || {}).label || "Value";
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
  sel.innerHTML = MAP_COLORS.map((m) => `<option value="${m.key}">${m.label}</option>`).join("");
  sel.value = state.mapColorBy;
  sel.onchange = () => { state.mapColorBy = sel.value; save(); if (state.view === "map") drawMap(); };

  const thr = $("mapThresh");
  thr.value = state.mapThreshold || 0;
  $("mapThreshVal").textContent = (+thr.value).toFixed(2);
  thr.oninput = () => { $("mapThreshVal").textContent = (+thr.value).toFixed(2); };
  thr.onchange = () => {
    state.mapThreshold = +thr.value;
    $("mapThreshVal").textContent = (+thr.value).toFixed(2);
    save(); if (state.view === "map") drawMap();
  };

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
        + `<div class="m">Composite <b>${isNum(r.__c) ? r.__c.toFixed(3) : "–"}</b>`
        + (state.mapColorBy !== "__c" ? ` · ${cl} <b>${isNum(v) ? v.toFixed(2) : "–"}</b>` : "") + `</div>`;
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
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.onclick = () => { loadOrg(b.dataset.org); });
  const xfer = $("xferToggle");
  if (xfer) {
    xfer.checked = state.orthoTransfer;
    xfer.onchange = () => {
      state.orthoTransfer = xfer.checked;
      applyTransferMode();
      save();
      recompute();
      if (_lastRow && $("drawer").classList.contains("open")) openDrawer(_lastRow);
    };
  }
  $("search").value = state.filters.search;
  let searchTimer = null;
  $("search").oninput = (e) => {
    state.filters.search = e.target.value.trim();
    if (searchTimer) clearTimeout(searchTimer);
    searchTimer = setTimeout(() => { save(); recompute(); }, 150);
  };
  $("resetWeights").onclick = () => {
    for (const c of COMPONENTS) state.weights[c.key] = { on: c.on, weight: c.weight };
    weightsDirty = true; buildWeights(); save(); recompute();
  };
  $("resetTiers").onclick = () => {
    for (const a of TIER_AXES) state.tiers[a.key] = [];
    buildTierPanel(); buildFilterBar(); save(); recompute();
  };
  $("colBtn").onclick = (e) => {
    e.stopPropagation();
    const open = $("colmenu").classList.toggle("open");
    if (open) {
      const pop = $("colPop"), r = $("colBtn").getBoundingClientRect();
      pop.style.top = Math.round(r.bottom + 4) + "px";
      const w = pop.offsetWidth || 210;
      const left = Math.max(8, Math.min(r.right - w, innerWidth - w - 10));
      pop.style.left = Math.round(left) + "px";
    }
  };
  document.addEventListener("click", (e) => { if (!$("colmenu").contains(e.target)) $("colmenu").classList.remove("open"); });
  $("exportBtn").onclick = exportCSV;
  $("methodsBtn").onclick = openMethods;
  $("methodsScrim").onclick = closeMethods;
  updateShortlistBtn();
  $("shortlistBtn").onclick = () => { state.shortlistOnly = !state.shortlistOnly; updateShortlistBtn(); save(); recompute(); };
  $("copyBtn").onclick = () => {
    const b = $("copyBtn");
    try { navigator.clipboard && navigator.clipboard.writeText(location.href); } catch (e) {}
    const t = b.textContent; b.textContent = "Copied ✓"; setTimeout(() => { b.textContent = t; }, 1200);
  };
  $("scrim").onclick = closeDrawer;
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") { closeDrawer(); closeMethods(); }
    if (state.view !== "map" && !e.target.matches("input,select,textarea")) {
      if (e.key === "ArrowRight") gotoPage(page + 1);
      if (e.key === "ArrowLeft") gotoPage(page - 1);
    }
  });
  addEventListener("resize", setTopbarOffset);
  loadOrg(state.org).then(() => {
    if (state._sel) { const r = DATA.rows.find((x) => x.uniprot_accession === state._sel); if (r) openDrawer(r); }
  });
}
init();
