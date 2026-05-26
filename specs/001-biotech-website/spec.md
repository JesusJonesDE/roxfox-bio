# Feature Specification: Biotech Company Website

**Feature Branch**: `001-biotech-website`

**Created**: 2026-05-26

**Status**: Draft

**Input**: User description: "A simple professional biotech company website for an AI-driven drug discovery startup with a pipeline of 3 drug candidates across two rare neurodegeneration disease areas."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - VC Investor Reviews the Company (Priority: P1)

A biotech venture capital investor receives the website URL before or after an introductory call. They need to quickly assess whether the company has a credible scientific foundation, a defined pipeline, and a clear value proposition worth pursuing further.

**Why this priority**: Investor credibility is the primary purpose of the site at this stage. Everything else is secondary.

**Independent Test**: A VC can land on the site, understand the company thesis, see the pipeline, and find an investor contact — all without leaving the page.

**Acceptance Scenarios**:

1. **Given** a VC visits the homepage, **When** they read the hero section, **Then** they immediately understand the company focuses on AI-driven rare neurodegeneration drug discovery
2. **Given** a VC scrolls to the pipeline section, **When** they view the table, **Then** they see 3 programs, their disease indications, modality, and pre-clinical stage — with no proprietary chemistry disclosed
3. **Given** a VC wants to reach out, **When** they reach the contact section, **Then** they find a working investor inquiry mechanism

---

### User Story 2 - Scientist or SAB Candidate Evaluates Credibility (Priority: P2)

A potential scientific advisor or academic collaborator is directed to the site to assess scientific legitimacy before agreeing to a meeting.

**Why this priority**: SAB recruitment is critical in the next 6 months and the website is a key signal.

**Independent Test**: A scientist with domain expertise reads the platform and pipeline sections and finds the framing credible and scientifically accurate.

**Acceptance Scenarios**:

1. **Given** a scientist visits the platform section, **When** they read it, **Then** they understand the AI discovery methodology at a high level without requiring proprietary implementation details
2. **Given** a scientist reviews the pipeline, **When** they see the disease indications, **Then** the framing is accurate and does not overclaim efficacy or clinical readiness

---

### User Story 3 - Founder Shares Site as a Credibility Signal (Priority: P3)

The founder shares the website URL in email introductions, grant applications, and professional profiles as proof that the company is operational.

**Why this priority**: A shareable, professional URL is a basic requirement for fundraising outreach.

**Independent Test**: The site loads correctly on mobile and desktop in under 3 seconds and looks professional enough to share without hesitation.

**Acceptance Scenarios**:

1. **Given** someone opens the URL on a mobile device, **When** the page loads, **Then** all content is readable and no elements are broken
2. **Given** the founder shares the link with a cold contact, **When** they view it, **Then** the design is consistent with established biotech company standards

---

### Edge Cases

- If the team section has only one person, it should display a single founder card with a note that the team is growing

- If a visitor searches for the program code names, no additional proprietary information should be discoverable beyond what is on the page

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Site MUST include a hero section with company name, one-line mission statement, and a primary call-to-action
- **FR-002**: Site MUST include a platform section describing the AI-driven discovery approach at a conceptual level — no proprietary methods or data
- **FR-003**: Site MUST include a pipeline table with columns: Program Code, Disease Indication, Modality, Stage
- **FR-004**: Pipeline table MUST use internal program codes (e.g. RXF-001) rather than target protein names, and MUST NOT display molecule structures, binding data, or any other proprietary chemistry
- **FR-005**: Site MUST include a team section; at minimum a single founder card plus a note that the team is being built
- **FR-006**: Site MUST include a contact/investor inquiry section with a direct email link (mailto) — no form required
- **FR-007**: Site MUST be a single page — no blog, no news, no pricing
- **FR-008**: Site MUST load in under 3 seconds on a standard broadband connection
- **FR-009**: Site MUST be fully responsive across desktop, tablet, and mobile viewports
- **FR-010**: Visual design MUST match the aesthetic standard of established AI biotech companies — clean, modern, scientific
- **FR-011**: Site MUST include complete SEO meta tags: `<title>`, `<meta name="description">`, Open Graph tags (`og:title`, `og:description`, `og:image`, `og:url`), and Twitter Card tags
- **FR-012**: Site MUST use semantic HTML5 structure (`<header>`, `<main>`, `<section>`, `<footer>`, proper heading hierarchy h1→h2→h3) to support search engine indexing
- **FR-013**: Site MUST include a `sitemap.xml` and `robots.txt` to guide crawlers

### Key Entities

- **Pipeline Program**: Internal code name, disease indication, modality (small molecule), current stage (Pre-Clinical)
- **Contact Inquiry**: Visitor name, email, inquiry type (Investor / Science / Other), message
- **Team Member**: Name, title, optional short bio

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A first-time visitor understands the company's disease focus within 10 seconds of landing
- **SC-002**: The full page loads in under 3 seconds on a standard broadband connection
- **SC-003**: All 3 pipeline programs are visible within the pipeline section without horizontal scrolling on desktop
- **SC-004**: The site renders correctly on screen widths from 375px to 1920px
- **SC-005**: An investor can locate the contact email within 60 seconds of arriving
- **SC-006**: Zero proprietary scientific data (structures, potency values, target names) is present anywhere on the site
- **SC-007**: Lighthouse SEO audit score ≥ 95
- **SC-008**: The company name and disease focus appear in the page `<title>` and meta description, making the site indexable for relevant search queries

## Assumptions

- Company name placeholder: **RoxFox Bio** — to be confirmed by founder before launch
- Program codes: **RXF-001** (SMA program), **RXF-002** (SMARD1 program), **RXF-003** (FTD program)
- Target protein names (VRK1, IGHMBP2, VCP) are NOT disclosed on the site until provisional patents are filed
- Disease indications are disclosed: Spinal Muscular Atrophy, SMARD1, Frontotemporal Dementia
- No backend required at launch — contact form uses a static form service or mailto link
- Hosting on a static platform (Vercel, Netlify, or GitHub Pages)
- No logo yet — typographic logo is acceptable for v1
- No analytics or cookie consent required at launch
