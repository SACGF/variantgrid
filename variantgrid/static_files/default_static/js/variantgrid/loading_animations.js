// @ts-check
// variantgrid/templates/default_templates/loading_animations.html
const stage = document.getElementById("stage");
const list = document.getElementById("loader-list");
const labelEl = document.getElementById("current-label");
// "white" = the analysis look: dark glyphs drawn over a clear (white) background.
// "dark"  = each loader paints its own dark background instead.
let stopCurrent = null, currentId = null, mode = "white";

function stop() {
  if (stopCurrent) { stopCurrent(); stopCurrent = null; }
  stage.innerHTML = "";
  currentId = null;
  labelEl.textContent = "-";
  [...list.children].forEach(b => b.removeAttribute("data-on"));
}

function play(id) {
  stop();
  const loader = VGLoaders.get(id);
  if (!loader) return;
  currentId = id;
  // Glyphs always use the "dark" palette; on the white background we let it show through
  // via clearBackground, exactly as the analysis overlay does.
  stopCurrent = VGLoaders.start(id, stage, {theme: "dark", clearBackground: mode === "white"});
  labelEl.textContent = loader.label;
  const btn = list.querySelector('[data-id="' + id + '"]');
  if (btn) btn.setAttribute("data-on", "");
}

VGLoaders.list().forEach(loader => {
  const btn = document.createElement("button");
  btn.textContent = loader.label;
  btn.dataset.id = loader.id;
  btn.addEventListener("click", () => play(loader.id));
  list.appendChild(btn);
});

function setMode(m) {
  mode = m;
  document.body.classList.toggle("dark", m === "dark");
  document.getElementById("theme-dark").toggleAttribute("data-on", m === "dark");
  document.getElementById("theme-light").toggleAttribute("data-on", m === "white");
  if (currentId) play(currentId);
}
document.getElementById("theme-light").addEventListener("click", () => setMode("white"));
document.getElementById("theme-dark").addEventListener("click", () => setMode("dark"));
document.getElementById("btn-random").addEventListener("click", () => play(VGLoaders.random().id));
document.getElementById("btn-stop").addEventListener("click", stop);

const first = VGLoaders.ids()[0];
if (first) play(first);
