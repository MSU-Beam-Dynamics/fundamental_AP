/* TrackPadWidgets — a tiny, dependency-free canvas plotting + slider engine.
 *
 * A widget is a set of pre-computed *frames*. Each frame is one point of a
 * (1- or 2-dimensional) parameter grid that Julia evaluated ahead of time with
 * TrackPad. Moving a slider only re-draws an already-computed frame, so the
 * published HTML needs no kernel.
 */
(function () {
  if (!window.__TPW__) {
    var TPW = {};
    window.__TPW__ = TPW;

    var STYLE = [
      '.tpw-root{font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;',
      'border:1px solid var(--tpw-line);border-radius:8px;padding:.75rem .9rem .9rem;margin:.4rem 0 1rem;',
      '--tpw-line:#d4d8de;--tpw-fg:#111827;--tpw-mid:#5b6472;--tpw-chip:#f3f4f6;}',
      '.tpw-root.tpw-dark{--tpw-line:#3f4652;--tpw-fg:#e5e7eb;--tpw-mid:#9aa4b2;--tpw-chip:#2b313b;}',
      '.tpw-title{font-weight:600;font-size:.95rem;color:var(--tpw-fg);margin-bottom:.5rem;}',
      '.tpw-ctls{display:flex;flex-wrap:wrap;gap:.4rem 1.2rem;align-items:center;margin-bottom:.6rem;}',
      '.tpw-ctl{display:flex;align-items:center;gap:.5rem;font-size:.82rem;color:var(--tpw-fg);}',
      '.tpw-ctl label{white-space:nowrap;}',
      '.tpw-ctl input[type=range]{width:170px;accent-color:#2563eb;cursor:pointer;}',
      '.tpw-val{font-variant-numeric:tabular-nums;min-width:5.2em;color:var(--tpw-mid);}',
      '.tpw-play{border:1px solid var(--tpw-line);background:var(--tpw-chip);color:var(--tpw-fg);',
      'border-radius:5px;padding:.1rem .5rem;font-size:.8rem;cursor:pointer;line-height:1.4;}',
      '.tpw-panels{display:flex;flex-wrap:wrap;gap:.6rem;}',
      '.tpw-panel{flex:1 1 300px;min-width:250px;position:relative;}',
      '.tpw-panel canvas{width:100%;display:block;}',
      '.tpw-outs{display:flex;flex-wrap:wrap;gap:.35rem .5rem;margin-top:.55rem;}',
      '.tpw-out{font-size:.78rem;background:var(--tpw-chip);color:var(--tpw-fg);',
      'border-radius:5px;padding:.12rem .45rem;font-variant-numeric:tabular-nums;}',
      '.tpw-out b{font-weight:600;color:var(--tpw-mid);font-variant-numeric:normal;}',
      '.tpw-note{font-size:.76rem;color:var(--tpw-mid);margin-top:.5rem;line-height:1.35;}'
    ].join('');
    var st = document.createElement('style');
    st.setAttribute('data-tpw', '1');
    st.textContent = STYLE;
    document.head.appendChild(st);

    /* ---------- helpers ---------------------------------------------- */

    function isDark() {
      var de = document.documentElement;
      if (de.classList.contains('dark')) return true;
      if (de.classList.contains('light')) return false;
      if (de.dataset && de.dataset.theme) return de.dataset.theme === 'dark';
      return !!(window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches);
    }

    function ticks(lo, hi, n) {
      n = n || 5;
      var span = hi - lo;
      if (!(span > 0)) return [lo];
      var step = Math.pow(10, Math.floor(Math.log(span / n) / Math.LN10));
      var err = span / n / step;
      if (err >= 7.5) step *= 10;
      else if (err >= 3.5) step *= 5;
      else if (err >= 1.5) step *= 2;
      var out = [], t = Math.ceil(lo / step - 1e-9) * step;
      for (; t <= hi + step * 1e-6; t += step) out.push(Math.abs(t) < step * 1e-6 ? 0 : t);
      out.step = step;
      return out;
    }

    function fmt(v, step) {
      if (v === 0) return '0';
      var a = Math.abs(v);
      if (a >= 1e5 || a < 1e-4) return v.toExponential(1);
      var d = Math.max(0, -Math.floor(Math.log(step || a) / Math.LN10));
      return v.toFixed(Math.min(d, 5));
    }

    /* --- logarithmic axes ------------------------------------------------
     * Positions on a log axis are carried in decades (log10 of the physical
     * value); only draw()'s data mapper converts, so callers keep passing
     * physical values.
     */
    var SUP = { '-': '\u207b', 0: '\u2070', 1: '\u00b9', 2: '\u00b2', 3: '\u00b3',
                4: '\u2074', 5: '\u2075', 6: '\u2076', 7: '\u2077', 8: '\u2078',
                9: '\u2079' };

    function sup(n) {
      return String(n).split('').map(function (c) { return SUP[c] || c; }).join('');
    }

    function logTicks(lo, hi) {
      var span = hi - lo, out = [];
      if (!(span > 0)) return [lo];
      if (span >= 2) {                       // whole decades, thinned if very wide
        var step = span > 14 ? 3 : (span > 7 ? 2 : 1);
        for (var d = Math.ceil(lo / step - 1e-9) * step; d <= hi + 1e-9; d += step) out.push(d);
      } else {                               // add 2 and 5 inside each decade
        var mant = [1, 2, 5];
        for (var e = Math.floor(lo); e <= Math.ceil(hi); e++) {
          for (var i = 0; i < 3; i++) {
            var v = e + Math.log(mant[i]) / Math.LN10;
            if (v >= lo - 1e-9 && v <= hi + 1e-9) out.push(v);
          }
        }
        out.sort(function (a, b) { return a - b; });
      }
      return out.length ? out : [lo, hi];
    }

    function logMinor(lo, hi) {
      var out = [];
      if (hi - lo > 9) return out;           // too dense to be useful
      for (var e = Math.floor(lo); e <= Math.ceil(hi); e++) {
        for (var m = 2; m <= 9; m++) {
          var v = e + Math.log(m) / Math.LN10;
          if (v >= lo && v <= hi) out.push(v);
        }
      }
      return out;
    }

    /* `plain` is decided once per axis so a single axis never mixes "1000" with
     * "10\u2076"; see plainLog() below. */
    function fmtLog(t, plain) {
      var k = Math.round(t), near = Math.abs(t - k) < 1e-6;
      if (plain) {
        var v = Math.pow(10, t);
        v = near ? Math.pow(10, k) : Math.round(v * 1e6) / 1e6;
        return v >= 1 ? String(v) : v.toFixed(Math.max(0, -Math.floor(t + 1e-9)));
      }
      if (near) return '10' + sup(k);
      var e = Math.floor(t + 1e-9);
      return Math.round(Math.pow(10, t - e)) + '\u00d710' + sup(e);
    }

    /* Plain decimal labels only while every tick stays in a readable range. */
    function plainLog(tk) {
      for (var i = 0; i < tk.length; i++) if (tk[i] < -3 || tk[i] > 3) return false;
      return true;
    }

    /* Limits from the series actually on screen, for panels that rescale per
     * frame. Mirrors _limits() in TrackPadWidgets.jl, including the 5% pad. */
    function dataExtent(series, key, log) {
      var lo = Infinity, hi = -Infinity;
      series.forEach(function (s) {
        var v = s[key];
        for (var i = 0; i < v.length; i++) {
          var t = v[i];
          if (t == null || !isFinite(t) || (log && !(t > 0))) continue;
          if (t < lo) lo = t;
          if (t > hi) hi = t;
        }
      });
      if (!isFinite(lo) || !isFinite(hi)) return log ? [1, 10] : [0, 1];
      if (log) {
        var l = Math.log(lo) / Math.LN10, h = Math.log(hi) / Math.LN10;
        var lp = (h - l < 1e-12) ? 0.5 : 0.05 * (h - l);
        return [Math.pow(10, l - lp), Math.pow(10, h + lp)];
      }
      if (hi - lo < 1e-12 * Math.max(1, Math.abs(hi))) {
        var q = Math.max(Math.abs(hi) * 0.05, 1);
        return [lo - q, hi + q];
      }
      var p = 0.05 * (hi - lo);
      return [lo - p, hi + p];
    }

    /* ---------- one panel --------------------------------------------- */

    /* `shared` carries the limits of the panel's share group, if it has one. */
    function draw(canvas, panel, series, dark, shared) {
      var css = getComputedStyle(canvas.parentNode);
      var fg = dark ? '#e5e7eb' : '#111827';
      var mid = dark ? '#9aa4b2' : '#5b6472';
      var grid = dark ? 'rgba(154,164,178,.20)' : 'rgba(91,100,114,.18)';
      var frame = dark ? '#5b6472' : '#9aa4b2';

      var dpr = window.devicePixelRatio || 1;
      var w = canvas.clientWidth || 300;
      var h = panel.height || 250;
      canvas.style.height = h + 'px';
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
      var g = canvas.getContext('2d');
      g.setTransform(dpr, 0, 0, dpr, 0, 0);
      g.clearRect(0, 0, w, h);
      g.font = '11px ' + (css.fontFamily || 'sans-serif');

      // Inset panels drop their numeric tick labels, and give the space back.
      var showtl = panel.ticklabels !== false;
      var hasY2 = panel.y2lim != null || series.some(function (s) { return s.axis === 'right'; });
      var m = { l: showtl ? 56 : (panel.ylabel ? 20 : 8),
                r: (hasY2 && showtl) ? 54 : 10,
                t: panel.title ? 20 : 8,
                b: showtl ? 34 : (panel.xlabel ? 18 : 8) };
      var pw = Math.max(10, w - m.l - m.r), ph = Math.max(10, h - m.t - m.b);

      var xlog = panel.xscale === 'log10', ylog = panel.yscale === 'log10';
      var L10 = function (v) {
        return (v == null || !(v > 0)) ? null : Math.log(v) / Math.LN10;
      };
      var TX = xlog ? L10 : function (v) { return v == null ? null : v; };
      var TY = ylog ? L10 : function (v) { return v == null ? null : v; };

      var xlim = panel.xlim || (shared && shared.x) || dataExtent(series, 'x', xlog);
      var ylim = panel.ylim || (shared && shared.y) || dataExtent(series, 'y', ylog);
      var x0 = TX(xlim[0]), x1 = TX(xlim[1]);
      var y0 = TY(ylim[0]), y1 = TY(ylim[1]);
      if (x0 == null || x1 == null) { x0 = 0; x1 = 1; }
      if (y0 == null || y1 == null) { y0 = 0; y1 = 1; }
      if (panel.equal && !xlog && !ylog) {
        // expand the shorter axis so one data unit is one pixel unit in both
        var sx = pw / (x1 - x0), sy = ph / (y1 - y0), s = Math.min(sx, sy);
        var cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
        x0 = cx - pw / (2 * s); x1 = cx + pw / (2 * s);
        y0 = cy - ph / (2 * s); y1 = cy + ph / (2 * s);
      }
      // Xr/Yr map already-transformed coordinates; X/Y take physical values.
      var Xr = function (t) { return m.l + (t - x0) / (x1 - x0) * pw; };
      var Yr = function (t) { return m.t + ph - (t - y0) / (y1 - y0) * ph; };
      var X = function (v) { var t = TX(v); return t == null ? null : Xr(t); };
      var Y = function (v) { var t = TY(v); return t == null ? null : Yr(t); };

      /* Right-hand axis: its own limits, its own linear mapping. */
      var y2lim = null, Y2 = Y;
      if (hasY2) {
        y2lim = panel.y2lim ||
                dataExtent(series.filter(function (s) { return s.axis === 'right'; }), 'y', false);
        var q0 = y2lim[0], q1 = y2lim[1];
        if (!(q1 > q0)) { q0 = q0 - 1; q1 = q1 + 1; }
        Y2 = function (v) {
          return (v == null || !isFinite(v)) ? null : m.t + ph - (v - q0) / (q1 - q0) * ph;
        };
      }

      // grid + ticks
      var tx = xlog ? logTicks(x0, x1) : ticks(x0, x1, Math.max(3, Math.round(pw / 90)));
      var ty = ylog ? logTicks(y0, y1) : ticks(y0, y1, Math.max(3, Math.round(ph / 55)));
      g.lineWidth = 1;
      if (xlog || ylog) {                     // faint 2…9 subdivisions inside each decade
        g.strokeStyle = grid;
        g.globalAlpha = 0.45;
        g.beginPath();
        if (xlog) logMinor(x0, x1).forEach(function (t) {
          var p = Math.round(Xr(t)) + .5; g.moveTo(p, m.t); g.lineTo(p, m.t + ph);
        });
        if (ylog) logMinor(y0, y1).forEach(function (t) {
          var p = Math.round(Yr(t)) + .5; g.moveTo(m.l, p); g.lineTo(m.l + pw, p);
        });
        g.stroke();
        g.globalAlpha = 1;
      }
      g.strokeStyle = grid;
      g.beginPath();
      tx.forEach(function (t) { var p = Math.round(Xr(t)) + .5; g.moveTo(p, m.t); g.lineTo(p, m.t + ph); });
      ty.forEach(function (t) { var p = Math.round(Yr(t)) + .5; g.moveTo(m.l, p); g.lineTo(m.l + pw, p); });
      g.stroke();

      g.strokeStyle = frame;
      g.strokeRect(m.l + .5, m.t + .5, pw, ph);

      if (showtl) {
        g.fillStyle = mid;
        g.textAlign = 'center'; g.textBaseline = 'top';
        var xplain = xlog && plainLog(tx), yplain = ylog && plainLog(ty);
        tx.forEach(function (t) {
          g.fillText(xlog ? fmtLog(t, xplain) : fmt(t, tx.step), Xr(t), m.t + ph + 5);
        });
        g.textAlign = 'right'; g.textBaseline = 'middle';
        ty.forEach(function (t) {
          g.fillText(ylog ? fmtLog(t, yplain) : fmt(t, ty.step), m.l - 6, Yr(t));
        });
      }

      if (hasY2 && showtl) {
        var t2 = ticks(y2lim[0], y2lim[1], Math.max(3, Math.round(ph / 55)));
        g.fillStyle = mid;
        g.textAlign = 'left'; g.textBaseline = 'middle';
        t2.forEach(function (t) { g.fillText(fmt(t, t2.step), m.l + pw + 6, Y2(t)); });
      }

      // axis labels
      g.fillStyle = fg;
      if (panel.xlabel) {
        g.textAlign = 'center'; g.textBaseline = 'bottom';
        g.fillText(panel.xlabel, m.l + pw / 2, h - 1);
      }
      if (panel.ylabel) {
        g.save();
        g.translate(11, m.t + ph / 2); g.rotate(-Math.PI / 2);
        g.textAlign = 'center'; g.textBaseline = 'top';
        g.fillText(panel.ylabel, 0, 0);
        g.restore();
      }
      if (panel.y2label) {
        g.save();
        g.translate(w - 2, m.t + ph / 2); g.rotate(Math.PI / 2);
        g.textAlign = 'center'; g.textBaseline = 'top';
        g.fillText(panel.y2label, 0, 0);
        g.restore();
      }
      if (panel.title) {
        g.textAlign = 'center'; g.textBaseline = 'top';
        g.fillText(panel.title, m.l + pw / 2, 2);
      }

      // data
      g.save();
      g.beginPath(); g.rect(m.l, m.t, pw, ph); g.clip();
      series.forEach(function (s) {
        var xs = s.x, ys = s.y, n = Math.min(xs.length, ys.length);
        var Yv = s.axis === 'right' ? Y2 : Y;
        g.strokeStyle = s.color; g.fillStyle = s.color;
        if (s.type === 'scatter') {
          var r = s.size || 2.6;
          g.globalAlpha = s.alpha == null ? 1 : s.alpha;
          for (var i = 0; i < n; i++) {
            var sx = X(xs[i]), sy = Yv(ys[i]);
            if (sx == null || sy == null) continue;
            g.beginPath(); g.arc(sx, sy, r, 0, 6.2832); g.fill();
          }
          g.globalAlpha = 1;
        } else {
          g.lineWidth = s.width || 1.8;
          var lineAlpha = s.alpha == null ? 1 : s.alpha;
          g.globalAlpha = lineAlpha;
          g.setLineDash(s.dash ? [5, 4] : []);
          /* Draw run by run, so a null breaks the stroke and any fill with it.
           * `fill` shades the run down to y = 0 before the outline is stroked. */
          var run = [];
          var flush = function () {
            if (run.length > 1) {
              if (s.fill) {
                var yz = Yv(0);
                if (yz != null) {
                  g.beginPath();
                  g.moveTo(run[0][0], yz);
                  for (var q = 0; q < run.length; q++) g.lineTo(run[q][0], run[q][1]);
                  g.lineTo(run[run.length - 1][0], yz);
                  g.closePath();
                  g.globalAlpha = lineAlpha * (s.fillAlpha == null ? 0.25 : s.fillAlpha);
                  g.fill();
                  g.globalAlpha = lineAlpha;
                }
              }
              g.beginPath();
              g.moveTo(run[0][0], run[0][1]);
              for (var r = 1; r < run.length; r++) g.lineTo(run[r][0], run[r][1]);
              g.stroke();
            }
            run = [];
          };
          for (var j = 0; j < n; j++) {
            var px = X(xs[j]), py = Yv(ys[j]);
            if (px == null || py == null) { flush(); continue; }
            run.push([px, py]);
          }
          flush();
          g.setLineDash([]);
          g.globalAlpha = 1;
        }
      });
      g.restore();

      // legend
      var leg = series.filter(function (s) { return s.label; });
      if (leg.length) {
        g.textAlign = 'left'; g.textBaseline = 'middle';
        var lw = 0;
        leg.forEach(function (s) { lw = Math.max(lw, g.measureText(s.label).width); });
        var bw = lw + 30, bh = leg.length * 14 + 8;
        var lp = panel.legend || 'right';
        var bx = /left/.test(lp) ? m.l + 6 : m.l + pw - bw - 6;
        var by = /bottom/.test(lp) ? m.t + ph - bh - 6 : m.t + 6;
        g.fillStyle = dark ? 'rgba(24,28,35,.78)' : 'rgba(255,255,255,.80)';
        g.strokeStyle = frame;
        g.beginPath(); g.rect(bx, by, bw, bh); g.fill(); g.stroke();
        leg.forEach(function (s, i) {
          var yy = by + 11 + i * 14;
          g.strokeStyle = s.color; g.fillStyle = s.color; g.lineWidth = 2;
          g.setLineDash(s.dash ? [4, 3] : []);
          if (s.type === 'scatter') { g.beginPath(); g.arc(bx + 12, yy, 3, 0, 6.2832); g.fill(); }
          else { g.beginPath(); g.moveTo(bx + 6, yy); g.lineTo(bx + 20, yy); g.stroke(); }
          g.setLineDash([]);
          g.fillStyle = fg;
          g.fillText(s.label, bx + 25, yy);
        });
      }
    }

    /* ---------- mount --------------------------------------------------- */

    TPW.mount = function (id, spec) {
      var root = document.getElementById(id);
      if (!root || root.dataset.tpwReady === '1') return;
      root.dataset.tpwReady = '1';
      root.className = 'tpw-root';

      var idx = spec.sliders.map(function (s) { return s.init || 0; });
      var canvases = [];

      if (spec.title) {
        var t = document.createElement('div');
        t.className = 'tpw-title'; t.textContent = spec.title;
        root.appendChild(t);
      }

      var ctls = document.createElement('div');
      ctls.className = 'tpw-ctls';
      root.appendChild(ctls);

      var panels = document.createElement('div');
      panels.className = 'tpw-panels';
      root.appendChild(panels);
      spec.panels.forEach(function (p) {
        var d = document.createElement('div');
        d.className = 'tpw-panel';
        if (p.basis) d.style.flexBasis = p.basis;
        if (p.minwidth != null) d.style.minWidth = p.minwidth + 'px';
        var c = document.createElement('canvas');
        d.appendChild(c); panels.appendChild(d); canvases.push(c);
      });

      var outs = document.createElement('div');
      outs.className = 'tpw-outs';
      root.appendChild(outs);

      if (spec.note) {
        var nd = document.createElement('div');
        nd.className = 'tpw-note';
        nd.textContent = spec.note;
        root.appendChild(nd);
      }

      function frame() {
        var k = 0, mul = 1;
        for (var i = 0; i < idx.length; i++) {
          k += idx[i] * mul;
          mul *= spec.sliders[i].values.length;
        }
        return spec.frames[k];
      }

      var vals = [];
      function render() {
        var dark = isDark();
        root.classList.toggle('tpw-dark', dark);
        var f = frame();
        spec.sliders.forEach(function (s, i) {
          vals[i].textContent = s.values[idx[i]];
        });
        var all = (spec.statics || []).concat(f.series || []);

        /* Panels naming the same share group get one set of limits, wide enough
         * for every point on any of them, recomputed for the frame on show. */
        var groups = {}, lims = {};
        spec.panels.forEach(function (p, i) {
          if (p.share) (groups[p.share] || (groups[p.share] = [])).push(i);
        });
        Object.keys(groups).forEach(function (name) {
          var members = groups[name];
          var ser = all.filter(function (s) { return members.indexOf(s.panel || 0) >= 0; });
          var p0 = spec.panels[members[0]];
          lims[name] = { x: dataExtent(ser, 'x', p0.xscale === 'log10'),
                         y: dataExtent(ser, 'y', p0.yscale === 'log10') };
        });

        spec.panels.forEach(function (p, i) {
          draw(canvases[i], p, all.filter(function (s) { return (s.panel || 0) === i; }),
               dark, p.share ? lims[p.share] : null);
        });
        outs.innerHTML = '';
        (f.readouts || []).forEach(function (r) {
          var e = document.createElement('span');
          e.className = 'tpw-out';
          e.innerHTML = '<b></b> ';
          e.firstChild.textContent = r[0] + ':';
          e.appendChild(document.createTextNode(r[1]));
          outs.appendChild(e);
        });
      }

      spec.sliders.forEach(function (s, i) {
        var d = document.createElement('div');
        d.className = 'tpw-ctl';
        var lab = document.createElement('label');
        lab.textContent = s.label;
        var inp = document.createElement('input');
        inp.type = 'range'; inp.min = 0; inp.max = s.values.length - 1;
        inp.step = 1; inp.value = idx[i];
        inp.setAttribute('aria-label', s.label);
        var v = document.createElement('span');
        v.className = 'tpw-val';
        vals.push(v);
        inp.addEventListener('input', function () { idx[i] = +inp.value; render(); });
        d.appendChild(lab); d.appendChild(inp); d.appendChild(v);
        ctls.appendChild(d);
      });

      if (spec.play !== false && spec.sliders.length) {
        var btn = document.createElement('button');
        btn.className = 'tpw-play';
        btn.type = 'button';
        btn.textContent = '▶ play';
        var timer = null, dir = 1;
        btn.addEventListener('click', function () {
          if (timer) { clearInterval(timer); timer = null; btn.textContent = '▶ play'; return; }
          btn.textContent = '❚❚ pause';
          var n = spec.sliders[0].values.length;
          var slider = ctls.querySelector('input');
          timer = setInterval(function () {
            idx[0] += dir;
            if (idx[0] >= n - 1) { idx[0] = n - 1; dir = -1; }
            else if (idx[0] <= 0) { idx[0] = 0; dir = 1; }
            slider.value = idx[0];
            render();
          }, spec.interval || 90);
        });
        ctls.appendChild(btn);
      }

      var ro = window.ResizeObserver ? new ResizeObserver(render) : null;
      if (ro) ro.observe(panels);
      new MutationObserver(render).observe(document.documentElement,
        { attributes: true, attributeFilter: ['class', 'data-theme'] });
      if (window.matchMedia) {
        var mq = window.matchMedia('(prefers-color-scheme: dark)');
        if (mq.addEventListener) mq.addEventListener('change', render);
      }
      render();
      // The canvas has no width until the browser has laid the page out.
      requestAnimationFrame(render);
    };
  }
})();
