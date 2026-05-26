# Data Model: Biotech Company Website

**Branch**: `001-biotech-website` | **Date**: 2026-05-26

*Note: This is a static site. There is no database. This document defines the content entities and their structure for use in the HTML template.*

## Entities

### Pipeline Program

Represents one drug discovery program displayed in the pipeline table.

| Field | Type | Value / Constraint |
|---|---|---|
| code | string | Internal code: RXF-001, RXF-002, RXF-003 |
| indication | string | Public disease name — NO target protein names |
| modality | string | "Small Molecule" for all three programs |
| stage | string | "Pre-Clinical" for all three programs |
| highlight | string | One-line scientific framing (optional) |

**Current values:**

| Code | Indication | Modality | Stage | Highlight |
|---|---|---|---|---|
| RXF-001 | Spinal Muscular Atrophy (SMA) | Small Molecule | Pre-Clinical | Kinase target with oncology expansion potential |
| RXF-002 | SMARD1 (SMA with Respiratory Distress) | Small Molecule | Pre-Clinical | Helicase modulator, highest genetic validation score |
| RXF-003 | Frontotemporal Dementia (FTD) | Small Molecule | Pre-Clinical | ATPase target, 140+ crystal structures available |

---

### Contact

Single email address displayed as a mailto link. No form, no backend.

| Field | Value |
|---|---|
| Display | `hello@roxfoxbio.com` (placeholder — update before launch) |
| Link | `mailto:hello@roxfoxbio.com` |

---

### Team Member

| Field | Type | Required | Notes |
|---|---|---|---|
| name | string | Yes | Full name |
| title | string | Yes | Role / title |
| bio | string | No | 1-2 sentences max |
| photo | image | No | Placeholder avatar if absent |

**Current values:**

| Name | Title | Bio |
|---|---|---|
| [Founder name] | Founder & CEO | Background in medical computer science and AI-driven drug discovery. |
| — | — | We are actively building our scientific team. |

## State Transitions

Not applicable — static site, no state management.

## Validation Rules

No form validation required — contact is a plain mailto link.
