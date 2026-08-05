import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import originalHandler from './site.js';

const BUNDLE_URL = 'https://raw.githubusercontent.com/byDenoso/AuditData/dfe7961f07e29358d0dfbecf3c3aa3f1e37c96d8/peer-research-atlas/preview/site-bundle.b64';
const BUNDLE_SHA256 = '184f27d1d6795a81deac330697b6eb725aa9e9b1de906580633feea74cb4d9bd';
let bundlePromise;

const TYPES = {
  '.html': 'text/html; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.svg': 'image/svg+xml; charset=utf-8',
};

function digest(value) {
  return createHash('sha256').update(value).digest('hex');
}

function cleanError(value) {
  return String(value?.message || value || 'falha desconhecida').replace(/[\r\n]+/g, ' ').slice(0, 360);
}

function securityHeaders(res) {
  res.setHeader('x-content-type-options', 'nosniff');
  res.setHeader('referrer-policy', 'strict-origin-when-cross-origin');
  res.setHeader('permissions-policy', 'camera=(), microphone=(), geolocation=()');
  res.setHeader('x-frame-options', 'SAMEORIGIN');
  res.setHeader('x-robots-tag', 'noindex, nofollow');
}

async function loadBundle() {
  bundlePromise ||= (async () => {
    const response = await fetch(BUNDLE_URL, { cache: 'force-cache' });
    if (!response.ok) throw new Error(`bundle HTTP ${response.status}`);
    const encoded = (await response.text()).trim();
    const compressed = Buffer.from(encoded, 'base64');
    const actual = digest(compressed);
    if (actual !== BUNDLE_SHA256) {
      throw new Error(`bundle SHA divergente: esperado ${BUNDLE_SHA256}, recebido ${actual}`);
    }
    const bundle = JSON.parse(gunzipSync(compressed).toString('utf8'));
    for (const path of ['index.html', 'favicon.svg', 'assets/styles.css', 'assets/app.js', 'assets/charts.js', 'assets/manifest-core.js', 'assets/snapshot-store.js', 'assets/view-model.js']) {
      if (typeof bundle[path] !== 'string') throw new Error(`bundle sem ${path}`);
    }
    return bundle;
  })();
  return bundlePromise;
}

async function serveStatic(res, requested) {
  const path = requested === '' || requested === '/' ? 'index.html' : requested.replace(/^\/+/, '');
  const allowed = /^(?:index\.html|favicon\.svg|assets\/(?:styles\.css|app\.js|charts\.js|manifest-core\.js|snapshot-store\.js|view-model\.js))$/;
  if (!allowed.test(path)) {
    res.statusCode = 404;
    securityHeaders(res);
    return res.end('Not found');
  }

  try {
    const bundle = await loadBundle();
    const extension = path.slice(path.lastIndexOf('.'));
    res.statusCode = 200;
    res.setHeader('content-type', TYPES[extension] || 'application/octet-stream');
    res.setHeader('cache-control', path === 'index.html' ? 'public, max-age=0, must-revalidate' : 'public, max-age=31536000, immutable');
    securityHeaders(res);
    res.end(bundle[path]);
  } catch (error) {
    bundlePromise = undefined;
    res.statusCode = 502;
    res.setHeader('content-type', 'text/plain; charset=utf-8');
    securityHeaders(res);
    res.end(`PEER preview bundle unavailable: ${cleanError(error)}`);
  }
}

export default async function handler(req, res) {
  const path = String(req.query?.path || '');
  if (path === '__manifest__') return originalHandler(req, res);
  return serveStatic(res, path);
}
