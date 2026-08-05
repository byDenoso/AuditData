import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';

const CHUNK_COMMIT = '7845259ca890483fa875a05809d1acbe32410da9';
const CHUNK_ROOT = `https://raw.githubusercontent.com/byDenoso/AuditData/${CHUNK_COMMIT}/peer-research-atlas/preview-v2`;
const CHUNK_URLS = Array.from({ length: 6 }, (_, index) => `${CHUNK_ROOT}/chunk-${index}.txt`);
const BUNDLE_SHA256 = '5e640e8a0c2adf5fd5a1e7e54b707b5366150fef2077f4effeef15c5268935b1';
const PRIMARY_MANIFEST = 'https://peer-studio-cosmology.vercel.app/data/studio-manifest.json';
const FALLBACK_MANIFEST = 'https://peer-studio-cosmology-kr7bclr2i-denosooo2-1701s-projects.vercel.app/data/studio-manifest.json';
const PRIMARY_EXPECTED_SHA = process.env.PEER_MANIFEST_EXPECTED_SHA256 || '';
const FALLBACK_EXPECTED_SHA = process.env.PEER_FALLBACK_MANIFEST_SHA256 || '';
let bundlePromise;
let lastValidSnapshot;

const REQUIRED_BUNDLE_PATHS = [
  'index.html',
  'favicon.svg',
  'assets/styles.css',
  'assets/app.js',
  'assets/charts.js',
  'assets/manifest-core.js',
  'assets/snapshot-store.js',
  'assets/view-model.js',
];

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
  return String(value?.message || value || 'falha desconhecida')
    .replace(/[\r\n]+/g, ' ')
    .slice(0, 360);
}

function securityHeaders(res) {
  res.setHeader('x-content-type-options', 'nosniff');
  res.setHeader('referrer-policy', 'strict-origin-when-cross-origin');
  res.setHeader('permissions-policy', 'camera=(), microphone=(), geolocation=()');
  res.setHeader('x-frame-options', 'SAMEORIGIN');
  res.setHeader('x-robots-tag', 'noindex, nofollow');
}

async function fetchChunk(url, index) {
  const response = await fetch(url, { cache: 'force-cache' });
  if (!response.ok) throw new Error(`bundle chunk ${index} HTTP ${response.status}`);
  const chunk = (await response.text()).trim();
  if (!chunk) throw new Error(`bundle chunk ${index} vazio`);
  return chunk;
}

async function loadBundle() {
  bundlePromise ||= (async () => {
    const chunks = await Promise.all(CHUNK_URLS.map(fetchChunk));
    const encoded = chunks.join('');
    const compressed = Buffer.from(encoded, 'base64');
    const actual = digest(compressed);
    if (actual !== BUNDLE_SHA256) {
      throw new Error(`bundle SHA divergente: esperado ${BUNDLE_SHA256}, recebido ${actual}`);
    }
    const decoded = gunzipSync(compressed).toString('utf8');
    const bundle = JSON.parse(decoded);
    for (const path of REQUIRED_BUNDLE_PATHS) {
      if (typeof bundle[path] !== 'string') throw new Error(`bundle sem ${path}`);
    }
    return bundle;
  })();
  return bundlePromise;
}

function inferSchema(manifest) {
  if (['peer-studio/manifest-v1', 'peer-studio/canonical-snapshot-v1'].includes(manifest?.schema)) {
    return manifest.schema;
  }
  if (
    !manifest?.schema &&
    manifest?.version &&
    manifest?.generated_at &&
    Array.isArray(manifest?.sources) &&
    Array.isArray(manifest?.series) &&
    manifest?.detection_status &&
    manifest?.cosmology
  ) {
    return 'peer-studio/canonical-snapshot-v1';
  }
  return null;
}

function validateManifest(manifest) {
  const errors = [];
  const schema = inferSchema(manifest);
  const sourceIds = new Set();
  const seriesIds = new Set();
  if (!schema) errors.push('schema canônico não reconhecido');
  if (!manifest?.version) errors.push('version ausente');
  if (!manifest?.generated_at || Number.isNaN(Date.parse(manifest.generated_at))) errors.push('generated_at inválido');
  if (!Array.isArray(manifest?.sources)) errors.push('sources inválido');
  if (!Array.isArray(manifest?.series)) errors.push('series inválido');
  if (!manifest?.detection_status || typeof manifest.detection_status !== 'object') errors.push('detection_status ausente');
  for (const source of manifest?.sources || []) {
    if (!source?.id || sourceIds.has(source.id)) errors.push('source id ausente ou duplicado');
    if (source?.id) sourceIds.add(source.id);
    if (source?.sha256 !== undefined && !/^[a-f0-9]{64}$/i.test(String(source.sha256))) {
      errors.push(`SHA-256 inválido em ${source?.id || 'source'}`);
    }
  }
  for (const series of manifest?.series || []) {
    if (!series?.id || seriesIds.has(series.id)) errors.push('series id ausente ou duplicado');
    if (series?.id) seriesIds.add(series.id);
    for (const sourceId of series?.source_ids || []) {
      if (!sourceIds.has(sourceId)) errors.push(`source_id quebrado: ${sourceId}`);
    }
    if (Array.isArray(series?.x) && Array.isArray(series?.y) && series.x.length !== series.y.length) {
      errors.push(`x/y inconsistente: ${series?.id || 'series'}`);
    }
  }
  return { valid: errors.length === 0, errors, schema };
}

async function readManifest(url, source, expectedSha) {
  const response = await fetch(url, { cache: 'no-store', headers: { accept: 'application/json' } });
  if (!response.ok) throw new Error(`${source} HTTP ${response.status}`);
  const raw = await response.text();
  const sha256 = digest(raw);
  if (expectedSha && sha256.toLowerCase() !== expectedSha.toLowerCase()) {
    throw new Error(`${source} SHA-256 divergente: ${sha256}`);
  }
  const manifest = JSON.parse(raw);
  const checked = validateManifest(manifest);
  if (!checked.valid) throw new Error(`${source} inválido: ${checked.errors.join('; ')}`);
  return { raw, sha256, schema: checked.schema, generatedAt: manifest.generated_at, source };
}

function emitManifest(res, snapshot, stale, reason = '') {
  res.statusCode = 200;
  res.setHeader('content-type', 'application/json; charset=utf-8');
  res.setHeader('cache-control', 'public, max-age=0, s-maxage=60, stale-while-revalidate=300');
  res.setHeader('x-peer-manifest-schema', snapshot.schema);
  res.setHeader('x-peer-manifest-sha256', snapshot.sha256);
  res.setHeader('x-peer-snapshot-source', snapshot.source);
  res.setHeader('x-peer-snapshot-stale', stale ? 'true' : 'false');
  res.setHeader('x-peer-snapshot-generated-at', snapshot.generatedAt);
  if (reason) res.setHeader('x-peer-snapshot-reason', cleanError(reason));
  res.setHeader(
    'access-control-expose-headers',
    'x-peer-manifest-schema, x-peer-manifest-sha256, x-peer-snapshot-source, x-peer-snapshot-stale, x-peer-snapshot-generated-at, x-peer-snapshot-reason',
  );
  securityHeaders(res);
  res.end(snapshot.raw);
}

async function serveManifest(res) {
  const failures = [];
  const candidates = [
    { url: PRIMARY_MANIFEST, source: 'primary', sha: PRIMARY_EXPECTED_SHA },
    { url: FALLBACK_MANIFEST, source: 'pinned-fallback', sha: FALLBACK_EXPECTED_SHA },
  ];
  for (const candidate of candidates) {
    try {
      const snapshot = await readManifest(candidate.url, candidate.source, candidate.sha);
      lastValidSnapshot = snapshot;
      return emitManifest(res, snapshot, candidate.source !== 'primary', failures[0]);
    } catch (error) {
      failures.push(error);
    }
  }
  if (lastValidSnapshot) return emitManifest(res, lastValidSnapshot, true, failures[0]);
  res.statusCode = 502;
  res.setHeader('content-type', 'application/json; charset=utf-8');
  res.setHeader('cache-control', 'no-store');
  securityHeaders(res);
  res.end(JSON.stringify({ error: 'PEER_MANIFEST_UNAVAILABLE', detail: failures.map(cleanError).join(' | ') }));
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
    res.setHeader(
      'cache-control',
      path === 'index.html' ? 'public, max-age=0, must-revalidate' : 'public, max-age=31536000, immutable',
    );
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
  if (path === '__manifest__') return serveManifest(res);
  return serveStatic(res, path);
}
