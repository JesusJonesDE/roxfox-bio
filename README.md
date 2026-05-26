# RoxFox Bio — Company Website

Static single-page website for RoxFox Bio, an AI-driven drug discovery startup.

## Quick Start

Open `site/index.html` in any browser — no build step needed.

## Deploy

See `specs/001-biotech-website/quickstart.md` for full Firebase setup and deployment instructions.

## Before Launch

Replace all placeholders:
- `[Founder Name]` in `site/index.html`
- `hello@roxfoxbio.com` with real company email
- `RoxFox Bio` with confirmed company name (if changed)
- `your-firebase-project-id` in `.firebaserc` and `.github/workflows/deploy.yml`
- Production domain in `site/sitemap.xml` and `site/robots.txt`

## Structure

```
site/
├── index.html     # Full single-page site
├── robots.txt     # SEO crawler instructions
├── sitemap.xml    # Sitemap for search engines
└── assets/        # Favicon, OG image

firebase.json      # Firebase Hosting config
.firebaserc        # Firebase project alias
.github/workflows/deploy.yml   # Auto-deploy on push to main
```
