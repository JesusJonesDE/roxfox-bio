# Implementation Plan: Biotech Company Website

**Branch**: `001-biotech-website` | **Date**: 2026-05-26 | **Spec**: [spec.md](spec.md)

**Input**: Feature specification from `/specs/001-biotech-website/spec.md`

## Summary

Build a single-page static biotech company website for an AI-driven drug discovery startup. The site must establish credibility with biotech VCs and scientists by presenting a 3-program pre-clinical pipeline (using internal program codes only — no target protein names until patents filed), an AI platform description, a founder team section, and an investor contact form. No backend. Deployed as a static site. Designed to the visual standard of Recursion / Relay Therapeutics.

## Technical Context

**Language/Version**: HTML5 + CSS3 + vanilla JavaScript (ES2022)

**Primary Dependencies**: Tailwind CSS v3 (via CDN — no build step), Google Fonts (Inter + DM Sans), Firebase Hosting, GitHub Actions

**Storage**: N/A — fully static, no database

**Testing**: Manual browser testing across Chrome, Safari, Firefox; Lighthouse audit for performance score ≥ 90

**Target Platform**: All modern browsers; responsive 375px–1920px

**Project Type**: Static single-page website

**Performance Goals**: Lighthouse performance score ≥ 90; page load < 3s on broadband

**Constraints**: No proprietary data on any page; no backend infrastructure; deploys automatically via GitHub Actions on push to `main`

**Scale/Scope**: Single HTML file + assets; ~5 sections; estimated 300–400 lines of HTML

## Constitution Check

*Constitution is a blank template (not yet filled in for this project). Applying default principles:*

- **Simplicity**: Chosen stack (HTML + Tailwind CDN) is the simplest approach that meets all requirements. ✓
- **No over-engineering**: No framework, no build pipeline, no backend. ✓
- **Deployable independently**: Single file can be previewed and deployed without tooling. ✓
- **IP safety**: Target names excluded from all content per spec FR-004. ✓

**GATE: PASSED** — no violations.

## Project Structure

### Documentation (this feature)

```text
specs/001-biotech-website/
├── plan.md          ← this file
├── research.md      ← Phase 0 output
├── data-model.md    ← Phase 1 output
├── quickstart.md    ← Phase 1 output
├── contracts/       ← Phase 1 output
└── tasks.md         ← Phase 2 output (from /speckit-tasks)
```

### Source Code (repository root)

```text
site/
├── index.html           # Full single-page site
├── sitemap.xml          # Single-URL sitemap for SEO crawlers
├── robots.txt           # Allow all crawlers
└── assets/
    ├── favicon.ico
    └── og-image.png     # 1200×630px Open Graph image for link sharing

firebase.json            # Firebase Hosting config (public dir: site/)
.firebaserc              # Firebase project alias
.github/
└── workflows/
    └── deploy.yml       # GitHub Actions: push to main → firebase deploy
README.md                # Setup and deploy instructions
```

**Structure Decision**: Single static file with inline Tailwind + CDN dependencies. All content in `site/index.html`. Assets folder for favicon and social sharing image only. No build step, no package.json, no node_modules.

## Complexity Tracking

No constitution violations — complexity tracking not required.
