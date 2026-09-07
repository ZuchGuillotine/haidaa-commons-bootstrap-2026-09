tes_screening_demo_inputs.csv - DEMO / PLACEHOLDER INPUT TABLE
==============================================================

WHAT THIS IS
This file exists ONLY so that tes_screening.py can be executed and its numbers
reported when the curated table tes_property_table.csv is not yet present in
this artifacts directory. It is a pipeline-exercise fixture, not a curated
property compilation. Every row carries "DEMO PLACEHOLDER" in its notes column.

DO NOT cite any number in this file as a property value for engineering or for
a literature claim. Use tes_property_table.csv (produced by the
property-data-normalizer role) once it exists. tes_screening.py prefers that
file automatically and reports which file it used.

HOW THE VALUES WERE CHOSEN
Central values are order-of-magnitude-correct figures of the kind reported in
the cited reviews and reports. They were entered from recollection of those
sources and were NOT re-checked against the source text row by row, so
individual digits may be wrong. The low/high columns are deliberately wide to
express that. The citation on each row is the origin of the general figure,
not a verified page-level locator.

Three disagreement pairs are encoded on purpose so that the source-selection
code path in tes_screening.py is exercised:
  1. Solar Salt thermal conductivity: 0.52 W/m-K (Sandia correlation) vs
     0.45 W/m-K (lower end of a review compilation). Literature scatter for
     molten salt conductivity is genuinely large; the two specific numbers here
     are illustrative.
  2. MgCl2-KCl-NaCl cost per tonne: 350 vs 220 USD/tonne. Indicative bulk
     chloride prices really do differ by roughly this factor across sources,
     mostly because assumed purity differs; the two numbers are illustrative.
  3. AlSi12 latent heat: 515 vs 560 kJ/kg. Reported enthalpies of fusion for
     near-eutectic Al-Si do scatter over roughly this range; neither number
     here is traced to a primary measurement.

SCHEMA
Identical to the curated table so that one parser serves both:
  material                  free text, groups rows
  family                    coarse family label
  property                  density | heat_capacity | thermal_conductivity |
                            melting_point | latent_heat |
                            min_operating_temperature |
                            max_operating_temperature | cost_per_tonne
                            (the parser also accepts common synonyms)
  value                     central value
  unit                      kg/m3, J/kg-K, W/m-K, C, kJ/kg, USD_per_tonne, ...
                            (the parser converts kJ/kg-K, g/cm3, kJ/kg,
                             USD_per_kg, K and a few others)
  temperature_C             temperature the value applies at; may be blank
  value_low, value_high     plausible interval; blank means "use the default
                            +/- percentage stated in the script"
  method_or_basis           how the number was obtained
  source_citation           author/year/venue
  source_locator            DOI, arXiv id, or report number
  notes                     free text; every row here says DEMO PLACEHOLDER
  disagreement_flag         "yes" when the row is one of several competing
                            rows for the same material+property

MULTIPLE ROWS FOR ONE MATERIAL+PROPERTY
Treated as competing sources. The Monte Carlo draws one of them uniformly per
draw, then samples inside that row's interval. If competing rows carry
different temperature_C values, the parser first keeps only the rows whose
temperature is nearest the material's own operating-window midpoint, and logs
the rows it dropped.

MATERIALS COVERED (8)
Solar Salt 60/40 NaNO3-KNO3; Hitec NaNO3-KNO3-NaNO2; MgCl2-KCl-NaCl eutectic;
Li2CO3-Na2CO3-K2CO3 eutectic; high-temperature castable concrete; silica sand
packed bed; dense alumina refractory; AlSi12 metallic PCM.

KNOWN LIMITATIONS OF THE FIXTURE
- Density for the packed bed is bulk density, for the salts it is melt density,
  and for concrete/alumina it is monolith density. Volumetric figures of merit
  are therefore not on a common system basis.
- min_operating_temperature is a system design choice, not a material property.
  It is carried here because the energy-density figure of merit needs a cold
  temperature; it should be treated as a scenario assumption.
- Property temperature dependence is ignored: one value per property is used
  across the whole swing.
- Costs are media costs only, in mixed and mostly unstated price-basis years.

SOURCES NAMED IN THE FILE
Zavoico 2001, Sandia National Laboratories, SAND2001-3835.
Gil et al. 2010, Renewable and Sustainable Energy Reviews 14(1),
  doi:10.1016/j.rser.2009.07.035.
Laing et al. 2012, Proceedings of the IEEE 100(2),
  doi:10.1109/JPROC.2011.2154290.
Vignarooban et al. 2015, Applied Energy 146,
  doi:10.1016/j.apenergy.2015.01.125.
Ho 2016, Applied Thermal Engineering 109,
  doi:10.1016/j.applthermaleng.2016.04.103.
Ding and Bauer 2021, Engineering 7(3), doi:10.1016/j.eng.2020.06.027.
Kenisarin 2010, Renewable and Sustainable Energy Reviews 14(3),
  doi:10.1016/j.rser.2009.11.011.
Mehos et al. 2017, NREL/TP-5500-67464.
Serrano-Lopez et al. 2013, Chemical Engineering and Processing,
  arXiv:1307.7343 - identifier not re-verified.
