/*
 * Illumina Read loader - an Illumina-style readout. 4 channel rows (A/C/G/T):
 * a fluorescence dot lights in the base's row, coloured by Illumina's 2-channel
 * dye identities (C red, T green, A red+green amber, G no dye / dim). A 5th row
 * is the base call, resolving N -> base to match the lit channel.
 * Adapted from claude/loading_animations/illumina-loader.html.
 */
(function () {
    "use strict";

    VGLoaders.register("illumina", {
        label: "Illumina Read",
        start: function (container, opts) {
            const SEQ = "ACGTTGCAAGTCCGATAGGCTTACAGGTCATGCCAATGGTCAGATCCGTTAACGGCATTGACGTACCTGAA";
            const CH = "ACGT", FS = 18;
            // Per-theme palette. dots = per-base dye colours (2-channel red/green scheme);
            // dark glows additively ('lighter'), light darkens the page ('multiply').
            const PAL = {
                dark: {
                    bg: "#07090c", dim: [120, 134, 148], brt: [222, 240, 236], blend: "lighter",
                    dots: { A: [240, 205, 70], C: [235, 80, 80], G: [86, 98, 112], T: [90, 205, 120] }
                },
                light: {
                    bg: "#eef1f5", dim: [176, 186, 196], brt: [24, 40, 54], blend: "multiply",
                    dots: { A: [200, 150, 20], C: [205, 45, 45], G: [120, 132, 144], T: [40, 160, 72] }
                }
            };
            const P = PAL[opts && opts.theme === "light" ? "light" : "dark"];
            const clearBg = !!(opts && opts.clearBackground);   // transparent canvas (shows page bg)
            const SPEED = 5;
            const T_N = 0.55, T_SETTLE = 1.6, FADE_IN = 0.16;
            const SLOW = 1;     // global playback rate (1 = original)
            const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
            const L = SEQ.length;

            const cv = document.createElement("canvas");
            cv.style.width = "100%";
            cv.style.height = "100%";
            cv.style.display = "block";
            container.appendChild(cv);
            const ctx = cv.getContext("2d");

            let cols, charW, leftX, rightX, chRowH, callY, chanY0, t0 = 0, glows = {}, gs, raf;

            const COL = [];
            for (let i = 0; i <= 16; i++) {
                const v = i / 16;
                COL[i] = `rgb(${P.dim.map((d, k) => Math.round(d + (P.brt[k] - d) * v)).join(",")})`;
            }
            const col = v => COL[v <= 0 ? 0 : v >= 1 ? 16 : Math.round(v * 16)];

            const baseOf = k => SEQ[(((k % L) + L) % L)];
            function makeGlow(rad, rgb) {
                const s = Math.ceil(rad * 2.4), g = document.createElement("canvas"); g.width = g.height = s;
                const c = g.getContext("2d"), grd = c.createRadialGradient(s / 2, s / 2, 0, s / 2, s / 2, s / 2);
                grd.addColorStop(0, `rgba(${rgb[0]},${rgb[1]},${rgb[2]},0.95)`);
                grd.addColorStop(0.4, `rgba(${rgb[0]},${rgb[1]},${rgb[2]},0.5)`);
                grd.addColorStop(1, `rgba(${rgb[0]},${rgb[1]},${rgb[2]},0)`);
                c.fillStyle = grd; c.fillRect(0, 0, s, s); return g;
            }
            function init() {
                const dpr = Math.min(2, devicePixelRatio || 1);
                const w = cv.clientWidth || container.clientWidth;
                const h = cv.clientHeight || container.clientHeight;
                cv.width = w * dpr; cv.height = h * dpr; ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
                ctx.font = `${FS}px ${getComputedStyle(cv).fontFamily}`;
                charW = ctx.measureText("G").width;
                leftX = 16 + 20; rightX = w - (16 + 20); cols = Math.floor((rightX - leftX) / charW);
                const callH = FS * 1.8, gap = 14;
                chRowH = Math.min(30, (h - 2 * 22 - callH - gap) / 4);
                chanY0 = (h - (4 * chRowH + gap + callH)) / 2 + chRowH / 2;
                callY = chanY0 - chRowH / 2 + 4 * chRowH + gap + callH / 2;
                const rad = Math.min(chRowH, charW * 1.6) * 0.5; gs = Math.ceil(rad * 2.4);
                for (const b of CH) glows[b] = makeGlow(rad, P.dots[b]);
                t0 = 0;
            }
            const chanY = r => chanY0 + r * chRowH;

            function frame(now) {
                raf = requestAnimationFrame(frame);
                if (!t0) t0 = now; const t = reduce ? 0 : (now - t0) / 1000 * SLOW;
                const sx = SPEED * charW * t, offset = sx % charW, baseIdx = Math.floor(sx / charW);

                if (clearBg) ctx.clearRect(0, 0, cv.width, cv.height);
                else { ctx.fillStyle = P.bg; ctx.fillRect(0, 0, cv.width, cv.height); }
                // Channel labels on both sides, styled the same as the called bases below.
                ctx.textBaseline = "middle"; ctx.fillStyle = col(0);
                ctx.font = `${FS}px ${getComputedStyle(cv).fontFamily}`;
                ctx.textAlign = "left";
                for (let r = 0; r < 4; r++) ctx.fillText(CH[r], leftX - 20, chanY(r));
                ctx.textAlign = "right";
                for (let r = 0; r < 4; r++) ctx.fillText(CH[r], rightX + 20, chanY(r));
                ctx.textAlign = "left";

                for (let j = -1; j <= cols + 1; j++) {
                    const k = baseIdx + j, x = leftX - offset + j * charW, age = (rightX - x) / (SPEED * charW);
                    let ch, vb;
                    // Hold an unresolved 'N', then resolve to the called base (GATC).
                    if (age < T_SETTLE) { ch = "N"; vb = age < T_N ? 1 : 1 - (age - T_N) / (T_SETTLE - T_N); }
                    else { ch = baseOf(k); vb = 0; }
                    ctx.fillStyle = col(vb); ctx.fillText(ch, x, callY);
                }
                ctx.globalCompositeOperation = P.blend;
                for (let j = -1; j <= cols + 1; j++) {
                    const k = baseIdx + j, x = leftX - offset + j * charW;
                    if (x < leftX - charW || x > rightX + charW) continue;
                    const age = (rightX - x) / (SPEED * charW), b = baseOf(k), r = CH.indexOf(b);
                    let a = Math.min(1, age / FADE_IN) * (0.55 + 0.45 * Math.exp(-age / 3));
                    if (a <= 0.02) continue;
                    ctx.globalAlpha = a;
                    ctx.drawImage(glows[b], x + charW / 2 - gs / 2, chanY(r) - gs / 2);
                }
                ctx.globalAlpha = 1; ctx.globalCompositeOperation = "source-over";
            }

            let resizeTimer;
            const onResize = () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(init, 150); };
            addEventListener("resize", onResize);

            init();
            raf = requestAnimationFrame(frame);

            return function stop() {
                cancelAnimationFrame(raf);
                removeEventListener("resize", onResize);
                clearTimeout(resizeTimer);
                cv.remove();
            };
        }
    });
})();
