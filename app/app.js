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
  tiers: {},              // axisKey -> "any" | "low" | "med" | "high"
  filters: { search: "", ranges: {}, cats: {}, bools: {}, families: [] },
  visibleCols: null,      // Set of column keys
  sort: { key: "__c", dir: -1 },
};

const cache = {};         // org -> parsed json
let DATA = null;          // current org payload
let AVAIL = {};           // component key -> bool (has data)
let filtered = [];        // current filtered+sorted rows
let page = 0;   // current page index (0-based)

// ---------- persistence ----------------------------------------------------
function save() {
  const s = {
    org: state.org, view: state.view, mapColorBy: state.mapColorBy,
    tiers: state.tiers, mapThreshold: state.mapThreshold,
    weights: state.weights,
    filters: { search: state.filters.search, ranges: state.filters.ranges,
      cats: Object.fromEntries(Object.entries(state.filters.cats).map(([k, v]) => [k, [...v]])),
      bools: state.filters.bools, families: state.filters.families },
    visibleCols: state.visibleCols ? [...state.visibleCols] : null,
    sort: state.sort,
  };
  try { localStorage.setItem(LS_KEY, JSON.stringify(s)); } catch (e) {}
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
    if (s.filters) {
      state.filters.search = s.filters.search || "";
      state.filters.ranges = s.filters.ranges || {};
      state.filters.cats = Object.fromEntries(Object.entries(s.filters.cats || {}).map(([k, v]) => [k, new Set(v)]));
      state.filters.bools = s.filters.bools || {};
      state.filters.families = Array.isArray(s.filters.families) ? s.filters.families : [];
    }
    state.visibleCols = s.visibleCols ? new Set(s.visibleCols) : null;
    state.sort = s.sort || state.sort;
  } catch (e) {}
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
       + (amt >= 58 ? "color:#fff;" : "color:var(--ink);");
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
  }
}

// ---------- filtering / sorting -------------------------------------------
function passFilters(row) {
  const f = state.filters;
  if (f.search) {
    const q = f.search.toLowerCase();
    const g = (row.gene || "").toLowerCase(), a = (row.uniprot_accession || "").toLowerCase();
    if (!g.includes(q) && !a.includes(q)) return false;
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
      return b && isNum(v) && v >= b.lo && v < b.hi;
    });
    if (!inAny) return false;
  }
  return true;
}
function recompute() {
  computeAll();
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
  if (isMap) { drawMap(); }
  else { renderThead(); renderPage(); }
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
  return ordered;
}
function tipAttr(text) {
  return text ? ` data-tip="${String(text).replace(/"/g, "&quot;")}"` : "";
}
const ROTATE_TYPES = new Set(["score", "pct", "plddt", "num", "int", "bool"]);
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
    + th("__c", "Composite", { tip: LEADING_COL_DESC.__c, rot: true, color: AXIS_COLORS.composite });
  for (const c of cols)
    h += th(c.key, c.label, { tip: c.desc, rot: true, color: colColor(c.key) });

  // section band: merge consecutive columns sharing a group
  const bandSeq = [{ key: "__c", group: "Composite" }, ...cols];
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
      + `<td class="colgene"${gtitle}><div class="gene">${row.name || row.gene || row.uniprot_accession}</div>`
      + `<div class="acc">${row.uniprot_accession}</div></td>`
      + `<td class="rc"><span class="heat cscore" style="${heatStyle(cc, AXIS_COLORS.composite)}">`
      + `${isNum(cc) ? cc.toFixed(2) : "–"}</span></td>`;
    for (const col of cols) h += `<td class="rc">${cellHTML(row, col)}</td>`;
    tr.innerHTML = h;
    tr.onclick = () => openDrawer(row);
    frag.appendChild(tr);
  }
  const tb = $("tbody");
  tb.innerHTML = "";
  tb.appendChild(frag);
  $("table").hidden = false;
  $("loading").hidden = true;
  window.scrollTo({ top: 0 });
  renderPager();
  setTopbarOffset();
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
function detailValueHTML(type, v) {
  if (type === "bool") return boolHTML(v);
  if (type === "class") return classBadgeHTML(v) + (FC_BY_ID[v] ? ` ${FC_BY_ID[v].label}` : "");
  if (type === "text") return (v === null || v === undefined || v === "") ? "–" : String(v).replace(/;/g, "; ");
  return fmt(type, v);
}
function openDrawer(row) {
  state._sel = row.uniprot_accession;
  document.querySelectorAll("#tbody tr").forEach((tr) =>
    tr.classList.toggle("sel", tr.dataset.acc === row.uniprot_accession));

  const comps = COMPONENTS.filter((c) => AVAIL[c.key] && state.weights[c.key] && state.weights[c.key].on
    && state.weights[c.key].weight > 0);
  const wsum = comps.reduce((s, c) => s + (isNum(row[c.key]) ? state.weights[c.key].weight : 0), 0);
  const contribHTML = comps.map((c) => {
    const v = row[c.key];
    const contrib = (isNum(v) && wsum > 0) ? (state.weights[c.key].weight * v / wsum) : null;
    const w = Math.round((contrib || 0) * 100);
    return `<div class="c" style="--cc:${colColor(c.key)}"><span class="cl">${c.label}</span>`
      + `<div class="cbarmini"><i style="width:${w}%"></i></div>`
      + `<span class="cvv">${isNum(v) ? v.toFixed(2) : "–"}</span></div>`;
  }).join("");

  const links = externalLinks(row).map((l) =>
    `<a href="${l.href}" target="_blank" rel="noopener">${l.label} ↗</a>`).join("");

  const groupColor = (title) => {
    const t = title.toLowerCase();
    if (t.startsWith("essential")) return AXIS_COLORS.essentiality;
    if (t.startsWith("ligand")) return AXIS_COLORS.ligandability;
    if (t.startsWith("structure")) return AXIS_COLORS.structure;
    if (t.startsWith("novelty")) return AXIS_COLORS.novelty;
    if (t.startsWith("cross")) return AXIS_COLORS.human_selective;
    return AXIS_COLORS.composite;
  };
  let groups = "";
  for (const g of DETAIL_GROUPS) {
    const rowsHTML = g.fields.filter(([k]) => DATA.columns.includes(k))
      .map(([k, lab, type]) => `<div class="drow"><span class="dk">${lab}</span>`
        + `<span class="dv">${detailValueHTML(type, row[k])}</span></div>`).join("");
    if (rowsHTML) groups += `<div class="dgroup" style="--dgc:${groupColor(g.title)}"><h3>${g.title}</h3>${rowsHTML}</div>`;
  }

  $("drawer").innerHTML = `
    <div class="dhead">
      <button class="close" id="drawerClose">×</button>
      <div class="g">${row.gene || row.uniprot_accession}</div>
      <div class="a">${row.uniprot_accession} · ${ORGANISM_META[state.org].name} ${ORGANISM_META[state.org].strain}</div>
      <div class="comp">
        <div class="cbar"><div class="track"><div class="fill" style="width:${((isNum(row.__c) ? row.__c : 0) * 100).toFixed(0)}%"></div></div>
        <span class="v">${isNum(row.__c) ? row.__c.toFixed(3) : "–"}</span></div>
        <div class="contrib">${contribHTML || '<span class="dk">No active components</span>'}</div>
      </div>
      <div class="links">${links}</div>
    </div>
    <div class="dbody">${groups}</div>`;
  $("drawerClose").onclick = closeDrawer;
  $("drawer").classList.add("open");
  $("drawer").setAttribute("aria-hidden", "false");
  $("scrim").classList.add("open");
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
  for (const r of DATA.rows) { const v = r[axisKey]; if (isNum(v) && v >= band.lo && v < band.hi) n++; }
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
      sw.onchange = () => { w.on = sw.checked; div.classList.toggle("off", !w.on); save(); recompute(); };
      rng.oninput = () => { w.weight = +rng.value; wv.textContent = rng.value; };
      rng.onchange = () => { save(); recompute(); };
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
  const header = ["rank", "uniprot_accession", "gene", "composite", ...cols.map((c) => c.key)];
  const lines = [header.join(",")];
  filtered.forEach((row, i) => {
    const vals = [i + 1, row.uniprot_accession, row.gene || "", isNum(row.__c) ? row.__c.toFixed(4) : ""];
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
async function loadOrg(org) {
  state.org = org;
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.setAttribute("aria-selected", b.dataset.org === org));
  $("table").hidden = true; $("loading").hidden = false; $("loading").textContent = "Loading…";
  if (!cache[org]) {
    const meta = ORGANISM_META[org];
    const res = await fetch(meta.file);
    if (!res.ok) { $("loading").textContent = `Failed to load ${meta.file} (${res.status})`; return; }
    cache[org] = await res.json();
  }
  DATA = cache[org];
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
  closeDrawer();
  recompute();
  save();
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
}
function renderMapLegend() {
  const lbl = (MAP_COLORS.find((m) => m.key === state.mapColorBy) || {}).label || "Value";
  const col = colColor(state.mapColorBy);
  $("maplegend").innerHTML = `<div class="lt">${lbl}</div>`
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
  wrap.addEventListener("mousemove", (e) => {
    const rect = wrap.getBoundingClientRect();
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
  initTooltips();
  initMap();
  initCollapsers();
  $("orgToggle").querySelectorAll("button").forEach((b) =>
    b.onclick = () => { loadOrg(b.dataset.org); });
  $("search").value = state.filters.search;
  $("search").oninput = (e) => { state.filters.search = e.target.value.trim(); save(); recompute(); };
  $("resetWeights").onclick = () => {
    for (const c of COMPONENTS) state.weights[c.key] = { on: c.on, weight: c.weight };
    buildWeights(); save(); recompute();
  };
  $("resetTiers").onclick = () => {
    for (const a of TIER_AXES) state.tiers[a.key] = [];
    buildTierPanel(); buildFilterBar(); save(); recompute();
  };
  $("colBtn").onclick = (e) => { e.stopPropagation(); $("colmenu").classList.toggle("open"); };
  document.addEventListener("click", (e) => { if (!$("colmenu").contains(e.target)) $("colmenu").classList.remove("open"); });
  $("exportBtn").onclick = exportCSV;
  $("scrim").onclick = closeDrawer;
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") closeDrawer();
    if (state.view !== "map" && !e.target.matches("input,select,textarea")) {
      if (e.key === "ArrowRight") gotoPage(page + 1);
      if (e.key === "ArrowLeft") gotoPage(page - 1);
    }
  });
  addEventListener("resize", setTopbarOffset);
  loadOrg(state.org);
}
init();
