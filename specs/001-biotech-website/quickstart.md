# Quickstart: Biotech Company Website

## Prerequisites

- A web browser (to preview locally)
- Node.js + Firebase CLI (`npm install -g firebase-tools`)
- A Firebase account (free Spark plan)
- A GitHub repository for this project

## Preview Locally

Open `site/index.html` directly in any browser. No server or build step required.

```bash
open site/index.html   # macOS
```

## First-Time Firebase Setup

```bash
firebase login
firebase init hosting
# → Public directory: site
# → Single-page app: No
# → Overwrite index.html: No
firebase deploy --only hosting
```

## GitHub Actions Auto-Deploy

Every push to `main` triggers a deploy automatically via `.github/workflows/deploy.yml`.

**One-time setup:**
1. Get a Firebase CI token: `firebase login:ci` → copy the token
2. In GitHub repo → Settings → Secrets → New secret
3. Name: `FIREBASE_TOKEN`, Value: paste the token

After that, `git push origin main` = live site update. No manual steps.

## Before Launch Checklist

- [ ] Replace `RoxFox Bio` with final company name (in index.html + meta tags)
- [ ] Replace contact email placeholder with real company email
- [ ] Add founder name to the team section
- [ ] Update `sitemap.xml` with the real domain URL
- [ ] Register domain and connect in Firebase Hosting console
- [ ] Run Lighthouse audit in Chrome DevTools → Performance ≥ 90, SEO ≥ 95

## Update Pipeline Table

When provisional patents are filed, update the pipeline table to add target protein names:

In `index.html`, find the pipeline section and add a "Target" column:
```html
<!-- Add after the Code column -->
<td>VRK1</td>   <!-- RXF-001 -->
<td>IGHMBP2</td> <!-- RXF-002 -->
<td>VCP</td>     <!-- RXF-003 -->
```
