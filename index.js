// Import stylesheets
import './style.css';

const w = 600;
const h = 600;

const canvas = document.createElement('canvas');

canvas.width = w;
canvas.height = h;

document.body.appendChild(canvas);

const ctx = canvas.getContext('2d');

const data = ctx.createImageData(w, h);

const points = [];

const pointCount = 100;

for (let i = 0; i < pointCount; i++) {
  const c = {
    r: Math.random(),
    g: Math.random(),
    b: Math.random(),
  };
  const point = {
    x: Math.random() * w,
    y: Math.random() * h,
    color: c,
  };
  points.push(point);
}

const colors = {};
const combine = (a, b) => Math.max(a, b);
const f = (x) => Math.max(0, Math.min(255, Math.floor(x * 255)));
const maxd = Math.sqrt(w * w + h * h);
for (let i = 0; i < pointCount; i++) {
  const point = points[i];
  for (let x = 0; x < w; x++) {
    for (let y = 0; y < h; y++) {
      const idx = 4 * (y * w + x);
      let r = 0;
      let g = 0;
      let b = 0;
      const dx = point.x - x;
      const dy = point.y - y;
      const d = Math.sqrt(dx * dx + dy * dy);
      const e = d / (maxd * pointCount);
      const e2 = 1000.0 / e;
      r += point.color.r * e2;
      g += point.color.g * e2;
      b += point.color.b * e2;
      const t = x + ',' + y;
      if (t in colors) {
        r = combine(r, colors[t].r);
        g = combine(g, colors[t].g);
        b = combine(b, colors[t].b);
      }
      colors[t] = {
        r,
        g,
        b,
      };
      data.data[idx + 0] = f(r);
      data.data[idx + 1] = f(g);
      data.data[idx + 2] = f(b);
      data.data[idx + 3] = f(1);
    }
  }
}

ctx.putImageData(data, 0, 0);
