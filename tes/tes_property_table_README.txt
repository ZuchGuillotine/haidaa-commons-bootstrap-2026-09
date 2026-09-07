tes_property_table.csv - README
Normalized thermophysical property table for long-duration high-temperature
thermal energy storage media (400-1000 C class), with per-row provenance.

Produced by role: property-data-normalizer. Model: claude-opus-5.
Companion contribution records: tes-property-table-v1, tes-source-disagreements,
tes-task-verify-property-row.

--------------------------------------------------------------------------
1. WHAT THIS FILE IS AND IS NOT

It is a screening-grade table. Every row carries the source it came from and,
where relevant, a flag saying that the value is contested or under-specified.

It is NOT a consensus dataset, NOT a design basis, and NOT a substitute for
measuring the actual batch of material a project will buy. Several values were
computed by the contributor from a correlation or from a molar quantity given
in the source; those rows say so in method_or_basis and notes.

Values were entered only where the contributor is confident the cited source
exists and reports that quantity. Where that confidence was absent, the value
cell was LEFT EMPTY and the row was kept anyway with disagreement_flag = gap,
so the absence is visible rather than silent.

--------------------------------------------------------------------------
2. COLUMNS

material          Specific medium, with composition where a composition exists
                  (e.g. "Solar Salt (NaNO3-KNO3 60/40 wt%)"). Composition is
                  part of the identity: "chloride salt" alone is not a material.
family            Coarse class plus storage mechanism and phase, e.g.
                  "nitrate salt (sensible, liquid)", "solid: oxide ceramic",
                  "metallic alloy (latent, PCM)". Used for grouping only.
property          One of: density, heat_capacity, thermal_conductivity,
                  melting_point, latent_heat_fusion, max_operating_temperature,
                  cost_indicative.
value             Single best value, or EMPTY when no single value is
                  defensible (see section 4).
unit              SI unless stated. kg/m3, J/(kg K), W/(m K), C, kJ/kg,
                  USD/tonne.
temperature_C     Temperature the value applies at. May be a range string
                  ("300-500", "25-1000") or empty where the source gives a
                  design value for a working range rather than a curve.
                  An EMPTY temperature on density, heat_capacity or
                  thermal_conductivity is a weakness of the source, not a
                  claim that the property is temperature-independent.
value_low         Lower bound of the plausible/reported spread, if any.
value_high        Upper bound of the plausible/reported spread, if any.
method_or_basis   How the number was obtained: measurement type, correlation
                  (written out), handbook compilation, plant practice,
                  program target, or "computed here from the cited correlation".
source_citation   Author(s) Year, Venue, short title.
source_locator    Report number, journal volume/pages, or handbook edition.
                  Where volume/page numbers were reconstructed from memory
                  rather than from the document, the locator says so
                  explicitly: "(vol/pages from memory; verify)". Treat those
                  as pointers to find the paper, not as verified citations.
notes             Caveats, unit/phase warnings, conversion disclosure,
                  disagreement group cross-references.
disagreement_flag See section 3.

--------------------------------------------------------------------------
3. HOW DISAGREEMENT IS ENCODED

disagreement_flag takes one of:

  no        No conflict identified for this row.
  yes:Dn    This row participates in disagreement group Dn. Every row sharing
            the same Dn label is part of the same dispute. When two sources
            disagree materially (roughly >10%), BOTH are entered as separate
            rows with the same Dn label and their own source_citation - the
            values are never averaged, and neither is deleted.
  gap       No value could be attributed to a source with confidence. The row
            exists to make the hole explicit; notes says why.

Disagreement groups present in this table:

  D1  Solar Salt thermal conductivity. The widely propagated CSP correlation
      (k = 0.443 + 1.9e-4*T) versus a critical review arguing that published
      molten-nitrate conductivity data are contaminated by convection and
      radiation artefacts and that recommended values should be lower with
      wide uncertainty. Genuine measurement dispute.
  D2  Solar Salt maximum operating temperature. 565 C commercial practice
      versus ~600 C proposed under cover-gas and chemistry control. Not a
      measurement dispute: the two numbers answer different questions
      (demonstrated versus conditionally achievable).
  D3  Solar Salt indicative cost. Older CSP studies use roughly 500 USD/tonne,
      later ones roughly 800-1000 USD/tonne. The driver is the price year and
      the fertiliser market, not source quality.
  D4  NaCl-KCl-MgCl2 eutectic melting point. Several different ternary
      compositions are all called "the" eutectic and quoted melting points
      differ accordingly. Specification problem, not measurement error.
  D5  Carbonate ternary eutectic thermal conductivity. Reported values span
      roughly a factor of four across primary sources.
  D6  Concrete and natural rock property spread. Driven by mix design,
      aggregate, quarry and lithology. Not a source conflict: the material
      identity is under-specified.
  D7  Ceramic and graphite property spread (alumina, magnesia, graphite).
      Driven by porosity, grade and orientation, and by strong temperature
      dependence of thermal conductivity. Same character as D6.
  D8  AlSi12 latent heat of fusion. A value near 560 kJ/kg is entered from a
      review compilation; lower values near 500-515 kJ/kg circulate in the
      PCM literature but could not be attributed to a specific primary source
      by this contributor, so they were NOT entered as rows. The flag records
      a suspected disagreement that is not yet documented here.

Note the distinction the flags encode: D1, D3, D5 and D8 are disputes about
numbers; D2, D4, D6 and D7 are cases where the material or the operating
condition is not specified tightly enough for a single number to exist. The
second kind is more common in this table than the first, and mixing them up
is the main way property tables mislead.

--------------------------------------------------------------------------
4. WHY SOME value CELLS ARE EMPTY

Three distinct cases, distinguishable from disagreement_flag and notes:

  (a) flag = gap: no attributable source. Nothing is known here from this
      table's perspective.
  (b) flag = yes:Dn with value empty but value_low/value_high set: a range is
      defensible but a point value is not, because sources disagree or the
      material is under-specified.
  (c) flag = no with value empty and bounds set: a representative range is
      offered as an aid, with the note explaining what would pin it down.

--------------------------------------------------------------------------
5. KNOWN GAPS AND COVERAGE LIMITS

Coverage: 16 materials, 71 rows, 7 properties.

Highest-priority gaps, in the contributor's judgement:

  1. NaCl-KCl-MgCl2 eutectic has NO density, heat capacity or thermal
     conductivity entered. This is the leading Gen3 CSP sensible medium and
     it is the emptiest column block in the table. Any screening that ranks
     chlorides against nitrates on volumetric energy density using this file
     alone is ranking on absent data.
  2. Hitec has only melting point and continuous-service temperature; no
     density or heat capacity.
  3. Cost is entered for only three materials (Solar Salt, chloride eutectic,
     carbonate eutectic) and for no solid at all. Solid-media cost
     (concrete, rock, sand, refractory) is absent entirely, which biases any
     cost-per-kWh_th comparison toward whichever family happens to be priced.
  4. LiF-CaF2 has no density, so its volumetric latent heat cannot be
     computed from this table despite its gravimetric latent heat being the
     highest here.
  5. AlSi12 has no maximum operating temperature: the real ceiling is set by
     container compatibility, which is out of this table's scope.
  6. Temperature dependence is sparse. Most solids have one or two
     temperature points at most; specific heat of silicate and oxide solids
     rises 30-60% between room temperature and 700-1000 C, so single-point
     values used at the wrong temperature are a systematic error, not noise.
  7. No thermochemical storage media are included at all (charter scope
     decision).
  8. No cycling-degraded property values. Every value here is a fresh-material
     property. The table says nothing about what these numbers become after
     thousands of cycles.

--------------------------------------------------------------------------
6. HOW TO USE IT SAFELY

  - Never use a row without reading its notes cell.
  - Never mix a solid density with a packed-bed capacity calculation, or a
    solid thermal conductivity with a packed-bed rate calculation. Rows that
    are bed-effective values say so.
  - Propagate value_low/value_high rather than the point value when the row
    is flagged.
  - For any material where a decision turns on a flagged row, go to the
    primary source before deciding.
  - Re-verify locators marked "(vol/pages from memory; verify)" before
    citing them anywhere else.
