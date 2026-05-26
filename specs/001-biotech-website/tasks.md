# Tasks: Biotech Company Website

**Input**: Design documents from `/specs/001-biotech-website/`

**Branch**: `001-biotech-website`

**Stack**: HTML5 + Tailwind CSS (CDN) + vanilla JS | Firebase Hosting | GitHub Actions

**Organization**: Tasks grouped by user story — each phase is independently deployable.

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Project scaffolding, Firebase config, CI/CD pipeline

- [x] T001 Create directory structure: `site/`, `site/assets/`, `.github/workflows/` at repo root
- [x] T002 Create `firebase.json` with public dir `site/`, rewrites for SPA disabled, cache headers for assets
- [x] T003 Create `.firebaserc` with project alias placeholder (`your-firebase-project-id`)
- [x] T004 [P] Create `.github/workflows/deploy.yml` — GitHub Actions workflow: on push to `main` → `firebase deploy --only hosting` using `FIREBASE_TOKEN` secret
- [ ] T005 [P] Create `site/assets/favicon.ico` — placeholder 32×32 favicon (dark navy with "R" lettermark)
- [ ] T006 [P] Create `site/assets/og-image.png` — 1200×630px Open Graph image: dark navy background, company name, tagline "AI-Driven Drug Discovery"
- [x] T007 Create `site/robots.txt` — allow all crawlers, point to sitemap
- [x] T008 Create `site/sitemap.xml` — single URL entry with placeholder domain `https://roxfoxbio.com`
- [ ] T009 Create `README.md` at repo root — setup instructions referencing `quickstart.md`

**Checkpoint**: Infrastructure ready. Firebase config present. CI/CD workflow wired. No HTML yet.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: The HTML skeleton and global styles that every section depends on

**⚠️ CRITICAL**: All user story phases build on this foundation

- [x] T010 Create `site/index.html` — full HTML5 document shell: `<!DOCTYPE html>`, correct lang, charset, viewport meta, Tailwind CDN link (`https://cdn.tailwindcss.com`), Google Fonts (Inter + DM Serif Display), placeholder `<title>` and meta description
- [x] T011 Add complete SEO meta block to `<head>` in `site/index.html`: `<title>`, `<meta name="description">`, canonical link, `og:title`, `og:description`, `og:image`, `og:url`, `og:type`, `twitter:card`, `twitter:title`, `twitter:description`, `twitter:image`
- [x] T012 Add JSON-LD Organization structured data `<script type="application/ld+json">` block to `<head>` in `site/index.html` — fields: name, url, description, foundingDate, sameAs
- [x] T013 Define Tailwind config inline `<script>` in `site/index.html` — custom colors: navy `#0A1628`, teal `#00D4C8`, slate `#64748B`; custom fonts: `display: 'DM Serif Display'`, `body: 'Inter'`
- [x] T014 Add `<header>` with sticky nav to `site/index.html` — logo (typographic "RoxFox Bio"), nav links anchoring to each section (Platform, Pipeline, Team, Contact), mobile hamburger toggle
- [x] T015 Add `<footer>` to `site/index.html` — company name, copyright year (JS-injected), tagline, email link
- [x] T016 Add smooth-scroll behavior and mobile nav toggle JavaScript to `site/index.html` (vanilla JS, inline `<script>` at bottom of body)

**Checkpoint**: HTML shell opens in browser. Nav and footer render. No content sections yet.

---

## Phase 3: User Story 1 — VC Investor Reviews the Company (Priority: P1) 🎯 MVP

**Goal**: A VC can land, understand the thesis, see the pipeline, and find contact — all on one page.

**Independent Test**: Open `site/index.html` in browser. Within 10 seconds identify disease focus. Scroll to pipeline table and confirm 3 programs visible. Find email link in contact section. No proprietary data visible anywhere.

### Implementation

- [x] T017 [US1] Add hero section to `site/index.html` — full-viewport dark navy background, `<h1>` with DM Serif Display ("Discovering Treatments for Rare Neurological Diseases"), one-line mission subheading, two CTAs: "Our Pipeline" (anchor scroll) and "Contact Us" (anchor scroll), subtle animated gradient or particle background via CSS
- [x] T018 [P] [US1] Add platform section to `site/index.html` — 3-column feature grid: "Genetic Evidence Scoring", "Structure-Based Screening", "ADMET Prediction"; each card with icon, heading, 2-sentence description; section heading "Our AI Platform"; no proprietary method names
- [x] T019 [US1] Add pipeline section to `site/index.html` — responsive table with columns: Program, Indication, Modality, Stage; rows for RXF-001 (SMA), RXF-002 (SMARD1), RXF-003 (FTD); all stages "Pre-Clinical" with a styled phase badge; note below table: "Target proteins disclosed upon patent filing"
- [x] T020 [US1] Add contact section to `site/index.html` — section heading "Get in Touch", two-line intro for investors and collaborators, prominent mailto link styled as a button (`hello@roxfoxbio.com`), secondary line for scientific inquiries

**Checkpoint**: US1 complete. Open in browser — VC journey fully functional. Lighthouse audit: Performance ≥ 90, SEO ≥ 95.

---

## Phase 4: User Story 2 — Scientist / SAB Candidate Evaluates Credibility (Priority: P2)

**Goal**: A scientist reads the platform and pipeline sections and judges them as credible and non-trivial.

**Independent Test**: Read platform section — methodology described at conceptual level, no overclaiming. Pipeline section uses correct disease names and accurate modality. No efficacy claims made.

### Implementation

- [x] T021 [US2] Expand platform section copy in `site/index.html` — add a 2-paragraph narrative below the 3-column grid: paragraph 1 describes the target identification approach (genetic evidence from population studies, structural biology databases); paragraph 2 describes the molecular screening approach (computational docking, ADMET filtering, validation roadmap)
- [x] T022 [P] [US2] Add "Why These Diseases" subsection to pipeline section in `site/index.html` — one sentence per program explaining the unmet need and why current treatments don't address it (SMA subtypes with normal SMN, SMARD1 zero treatments, FTD no approved disease-modifying therapy)
- [x] T023 [P] [US2] Add scientific disclaimer / pipeline stage legend below pipeline table in `site/index.html` — define Pre-Clinical, clarify that programs are computational discovery stage, include IND timeline language ("targeting IND-enabling studies 2026–2027")

**Checkpoint**: US2 complete. A scientist reading the page finds accurate, credible scientific framing.

---

## Phase 5: User Story 3 — Founder Shares Site as Credibility Signal (Priority: P3)

**Goal**: Site renders perfectly on mobile, loads fast, looks professional enough to share cold.

**Independent Test**: Open on iPhone SE (375px width) — all text readable, no overflow. Run Lighthouse: Performance ≥ 90. Share URL on LinkedIn — OG image and title appear correctly in preview.

### Implementation

- [x] T024 [US3] Add team section to `site/index.html` — single founder card (name placeholder, title "Founder & CEO", 2-line bio), plus a "We're building the team" card with "Open Positions" placeholder text; section heading "Our Team"
- [x] T025 [P] [US3] Audit and fix all responsive breakpoints in `site/index.html` — test hero, platform grid, pipeline table, team cards at 375px, 768px, 1280px; pipeline table converts to card layout on mobile
- [x] T026 [P] [US3] Add `<meta name="theme-color">` and apple touch icon references to `<head>` in `site/index.html` for mobile browser chrome theming
- [x] T027 [P] [US3] Verify all Open Graph and Twitter Card meta tags render correctly — test with a meta tag validator (ogp.me or similar); confirm og:image resolves to absolute URL

**Checkpoint**: US3 complete. Site is shareable, mobile-ready, and passes social preview tests.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Performance, accessibility, final SEO, deploy verification

- [ ] T028 [P] Replace placeholder company name `RoxFox Bio` with confirmed final name throughout `site/index.html` (appears in nav, hero, footer, meta tags — use find/replace)
- [ ] T029 [P] Replace `hello@roxfoxbio.com` placeholder with real company email in `site/index.html` and `site/sitemap.xml`
- [ ] T030 [P] Update `site/sitemap.xml` with real production domain once registered
- [ ] T031 Add `loading="lazy"` to any images; minify inline CSS if any custom styles added; verify Tailwind CDN is loaded from correct version pin
- [ ] T032 Run Lighthouse audit in Chrome DevTools against local file — record scores for Performance, SEO, Accessibility, Best Practices; fix any issues scoring below 90/95
- [ ] T033 Validate HTML via W3C validator — fix any errors; confirm heading hierarchy is h1→h2→h3 throughout
- [ ] T034 Test full deploy pipeline: commit to `main` branch → confirm GitHub Actions workflow triggers → confirm Firebase deploy succeeds → confirm live URL serves correct content
- [ ] T035 [P] Add founder name and real bio to team section in `site/index.html`
- [ ] T036 Verify `robots.txt` and `sitemap.xml` are accessible at production domain root after deploy

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundation)**: Depends on Phase 1
- **Phase 3 (US1)**: Depends on Phase 2 — **MVP stops here**
- **Phase 4 (US2)**: Depends on Phase 2 — can run in parallel with Phase 3
- **Phase 5 (US3)**: Depends on Phase 2 — can run in parallel with Phase 3+4
- **Phase 6 (Polish)**: Depends on all content phases complete

### User Story Dependencies

- **US1**: Independent after Foundation — delivers the full VC journey
- **US2**: Independent after Foundation — expands scientific credibility copy
- **US3**: Independent after Foundation — mobile polish + team section

### Parallel Opportunities

```
# Phase 1 — run together:
T004 (GitHub Actions)  +  T005 (favicon)  +  T006 (OG image)  +  T007 (robots.txt)

# Phase 3 — run together after T017 (hero):
T018 (platform section)  +  T019 (pipeline section)  +  T020 (contact section)

# Phase 4 — run together after T021:
T022 (disease rationale)  +  T023 (stage legend)

# Phase 6 — run together:
T028 (name)  +  T029 (email)  +  T030 (sitemap)  +  T035 (team bio)
```

---

## Implementation Strategy

### MVP (Phase 1 + 2 + 3 only)

1. Complete Phase 1: Scaffolding + Firebase config
2. Complete Phase 2: HTML shell + nav + SEO head
3. Complete Phase 3: Hero + Platform + Pipeline + Contact
4. **STOP and VALIDATE**: Open in browser, run Lighthouse, confirm VC journey works
5. Deploy to Firebase: `firebase deploy --only hosting`
6. Share URL

**Time estimate**: 3–4 hours for MVP.

### Full Site

Add Phase 4 (scientific copy) + Phase 5 (mobile + team) + Phase 6 (polish) for a production-ready site.

**Total time estimate**: 6–8 hours.

---

## Notes

- No tests required — static marketing site, validated by browser + Lighthouse
- All placeholder values (company name, email, Firebase project ID) collected in `quickstart.md` pre-launch checklist
- Target names (VRK1, IGHMBP2, VCP) intentionally absent — add after provisional patents filed per `quickstart.md`
- Tailwind CDN is acceptable for this use case (single page, low traffic) — migrate to Tailwind CLI build if site grows
