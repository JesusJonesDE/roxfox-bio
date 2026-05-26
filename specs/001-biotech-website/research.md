# Research: Biotech Company Website

**Branch**: `001-biotech-website` | **Date**: 2026-05-26

## Decision Log

### D-001: Tech Stack

**Decision**: HTML5 + Tailwind CSS v3 (CDN) + vanilla JavaScript

**Rationale**: The site is a static marketing page with no dynamic data requirements. A build pipeline (Next.js, Astro, Vite) adds deployment complexity and maintenance overhead without benefit. Tailwind via CDN allows professional design without a `npm install`. The founder can edit the single HTML file directly.

**Alternatives considered**:
- Next.js — overkill, requires Node.js deployment, adds complexity
- Astro — good option but requires build step; founder may not have Node.js tooling
- Pure CSS — slower to produce professional design; Tailwind utilities are faster

---

### D-002: Contact Mechanism

**Decision**: Plain `mailto:` link — no form

**Rationale**: Simplest possible approach. No backend, no third-party service, no maintenance. Investors and scientists who want to reach out will use email directly. A mailto link is sufficient for a pre-launch credibility site.

---

### D-003: Hosting

**Decision**: Firebase Hosting + GitHub Actions CI/CD

**Rationale**: Firebase Hosting provides a global CDN, automatic HTTPS, custom domain support, and fast deploy times. GitHub Actions automates deployment on every push to `main` — no manual steps after initial setup. Free Spark plan covers static hosting at this scale indefinitely.

**Setup required**:
1. Create Firebase project in Firebase Console
2. Run `firebase init hosting` locally to generate `firebase.json` and `.firebaserc`
3. Add `FIREBASE_TOKEN` secret to GitHub repository settings
4. GitHub Actions workflow: on push to `main` → `firebase deploy --only hosting`

**Alternatives considered**:
- Netlify — good but Firebase is preferred by the founder
- Vercel — excellent for Next.js but no advantage for plain HTML
- GitHub Pages — no CDN, no custom redirect rules

---

### D-004: Typography

**Decision**: Inter (body + UI) + DM Serif Display (hero headline only)

**Rationale**: Inter is the de facto standard for clean, scientific tech/biotech UIs (used by Recursion, Relay, Exscientia). DM Serif adds one premium typographic moment in the hero without feeling pharmaceutical-generic.

---

### D-005: Color Palette

**Decision**: Deep navy (#0A1628) background sections + white (#FFFFFF) primary + electric teal accent (#00D4C8) + slate gray (#64748B) secondary text

**Rationale**: Deep navy reads as serious/scientific and is used by top-tier biotech companies (Blueprint Medicines, Relay Therapeutics). Teal accent avoids the cliché blue-on-blue pharma look and nods to data/AI context. High contrast for accessibility.

---

### D-006: Pipeline Disclosure Level

**Decision**: Disease indications disclosed; target protein names withheld; program codes used (RXF-001, RXF-002, RXF-003)

**Rationale**: Provisional patents not yet filed. Naming VRK1, IGHMBP2, VCP publicly before patent filing could compromise novelty. Disease names (SMA, SMARD1, FTD) are public knowledge and provide sufficient credibility signal without IP risk. This follows standard early-stage biotech practice.

---

### D-007: Company Name

**Decision**: **RoxFox Bio** (placeholder from session name)

**Rationale**: Used as working title throughout build. Founder to confirm before domain registration and launch. Easy to replace — appears in 3 places in the HTML.

---

### D-009: SEO Strategy

**Decision**: Full on-page SEO — semantic HTML, meta tags, Open Graph, sitemap.xml, robots.txt, structured data (Organization schema)

**Rationale**: A static HTML site is inherently SEO-friendly if built correctly. Key targets: searches for the company name, founder name, "AI drug discovery SMA", "VRK1 inhibitor" (after patent filing). Open Graph tags ensure the site looks professional when shared on LinkedIn by investors/founders.

**What to include**:
- `<title>`: "RoxFox Bio — AI-Driven Drug Discovery for Rare Neurodegeneration"
- `<meta name="description">`: 155-char summary covering AI platform + SMA/FTD
- `og:image`: 1200×630px branded image (generated as `assets/og-image.png`)
- JSON-LD Organization schema: name, URL, description, foundingDate
- `sitemap.xml`: single URL entry
- `robots.txt`: allow all crawlers

---

### D-008: AI Platform Description Copy

**Decision**: Describe methodology as "structure-based virtual screening + ADMET prediction + genetic evidence scoring" without naming specific tools (AlphaFold, DiffDock, Open Targets)

**Rationale**: The tools are open-source and not IP-sensitive, but naming them in marketing copy sounds academic rather than proprietary. Abstract description ("AI-driven target identification and molecular screening platform") is stronger for investor communication.
