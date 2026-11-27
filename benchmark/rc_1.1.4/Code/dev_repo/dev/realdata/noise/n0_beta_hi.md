# N0: physiological upper limit β_hi of the missegregation rate (literature)

Written 2026-11-20, following N0 of `PLAN.md`. All units below are "per chromosome copy per division" (KARST's β).

**Conversion:**
- Imaging gives "the fraction of divisions with an error"; assuming each error involves only one chromosome, divide by the number of copies per cell (about 46).
- MDD (missegregations per diploid division) divided by 46.

**Verification:**
- Every reference was checked by a sub-agent on Europe PMC / PubMed; for those marked FT, the numbers were read in the full text.
- I also spot-checked two of them, Thompson & Compton 2008 and Lynch et al. 2024; titles, authors, journals and years are all correct.

## Recommendation: β_hi = 1×10⁻² (200 times the default 5e-5)

That is, about 0.5 missegregations per division with 46 copies. Basis:
- **CIN cancer cell lines:** Caco2 and MCF-7 keep proliferating at missegregation rates of about 0.7–1% per chromosome per division (Thompson & Compton 2008, FISH); U2OS about 1% (Bakhoum et al. 2009; Lynch et al. 2024, 0.33–0.46 MDD).
- **Tolerance experiments:** after knockdown of Mps1/BubR1 in tumour cell lines, more than 50% of divisions are erroneous and colony formation is still unaffected; when the errors are further increased with paclitaxel, colonies can no longer form (Janssen et al. 2009).
- **Patient tumours:** in 54 DLBCL biopsies, 39.5 ± 18.1% of anaphases have lagging chromosomes or chromosome bridges, about 8.6×10⁻³. The authors call "one event every two divisions" the high end of the tolerable range (Bakhoum et al. 2011).
- **Lethality:** missegregating 5 or more chromosomes per division is usually lethal, about 0.1 (Funk et al. 2016 review; Silk et al. 2013 in mouse cells: with 34% of divisions erroneous, mostly single chromosomes, it is tolerated; with 50%, mostly multiple chromosomes, the cells die).

## Values by source

| source | system, method | raw value | β (per copy per division) |
|---|---|---|---|
| Thompson & Compton 2008, J Cell Biol (FT) | RPE1, HCT116; FISH on daughter-cell pairs | about 0.025% per chromosome | 2.5×10⁻⁴ |
| Lynch et al. 2024, PNAS (FT) | RPE1 | 0.01–0.05 MDD | 2×10⁻⁴ – 1.1×10⁻³ |
| Bolhaqueiro et al. 2019, Nat Genet (FT) | normal colon organoids; live-cell imaging | 1.5–6.1% of divisions erroneous | 3×10⁻⁴ – 1.3×10⁻³ |
| Lengauer et al. 1997, Nature | colorectal cancer cell lines; clonal FISH | > 10⁻² per chromosome per generation | ≳ 10⁻² |
| Thompson & Compton 2008 (FT) | HT29, MCF-7, Caco2; FISH | 0.3–1.0% per chromosome | 3×10⁻³ – 1×10⁻² |
| Bakhoum et al. 2009, Nat Cell Biol (FT) | U2OS; FISH | about 1% per chromosome | about 1×10⁻² |
| Lynch et al. 2024 (FT) | U2OS | 0.33–0.46 MDD | 7×10⁻³ – 1×10⁻² |
| Bolhaqueiro et al. 2019 (FT) | 11 colorectal tumour organoids; imaging | about 10–50% of divisions erroneous | 2×10⁻³ – 1.1×10⁻² |
| Bollen et al. 2021, Nat Genet (FT) | colorectal tumour organoids; imaging, single-cell sequencing, ABC | about 0.1 whole-chromosome CNA per division | about 2×10⁻³ |
| Lynch et al. 2022, eLife (FT) | organoids, breast tumour scDNA; ABM + ABC | Pmisseg | organoids 9×10⁻⁴ – 2.2×10⁻²; breast tumours 5×10⁻³ – 3×10⁻² |
| Bakhoum et al. 2011, Clin Cancer Res (FT) | 54 DLBCL biopsies; anaphase | 39.5 ± 18.1% abnormal anaphases | about 8.6×10⁻³ |
| Janssen et al. 2009, PNAS (FT) | tumour cell lines, Mps1/BubR1 knockdown | > 50% of divisions erroneous, still proliferating | tolerable ≳ 1×10⁻² |
| Silk et al. 2013, PNAS (FT, mouse) | mouse fibroblasts | 34% tolerated; 50% (multiple chromosomes) lethal | about 9×10⁻³ vs ≳ 2×10⁻² |
| Komarova & Wodarz 2004; Elizalde et al. 2018 | theory | optimum 10⁻³–10⁻²; survival peaks below 10⁻³ | 10⁻³ – 10⁻² |

- **No usable numbers:** Salehi et al. 2021 (no missegregation rate in the full text); Funnell et al. 2022, Kimmel et al. 2023 (only in figures); CINner (Dinh et al. 2025) takes β = 5×10⁻⁵ as an assumption, not a measurement.
- **Original text could not be verified:** the numbers of Laughney et al. 2015 (about 1.3×10⁻³, theoretical optimum about 1.9×10⁻³), seen only in secondary reviews.
- **PDX:** no publication gives a per-copy missegregation rate for serially passaged PDX. Current evidence shows that raw missegregation rates in vivo and in vitro are of the same order of magnitude; tumour karyotype diversity is lower than imaging would predict, which Lynch et al. 2022 attribute to selection.

## Uncertainty

- **The conversion has an error of 2–4 times:**
  - Imaging overestimates: most lagging chromosomes eventually still segregate correctly (Thompson & Compton 2011), and chromosome bridges are mostly structural.
  - FISH underestimates: it misses micronuclei; the rate from one probe may correspond to 2–4 copies.
  - Model-based ABC estimates vary 2–5 times with the growth model.
- **These are all raw missegregation rates, with selection not subtracted.** For lineages that survive, the effective rate is lower. For the purpose of an "upper limit", the raw rate is a suitable upper bound.

## Implications for this study (for reference in your decision)

1. **KARST's default 5×10⁻⁵ is much lower than the measured values:** about 5 times lower than the stable diploid RPE1, and 20–200 times lower than CIN cell lines. But the β study showed that raising the β used in the fit brings no benefit, either in simulations or on Salehi. So the default is more like an "effective" parameter than a physiological quantity.
2. **The best forecast missegregation rate in the decomposition study lies within the measured range.** When only the forecast missegregation rate is changed, the best is ×16–×128, i.e. 8×10⁻⁴ – 6.4×10⁻³, which falls exactly in the measured range of CIN cell lines and tumour organoids.
3. **With β_hi = 1×10⁻², the noise flag is hard to trigger.** r̂ must exceed 200 times the default to count as "noise-dominated". Because r̂ is inversely proportional to the assumed division rate, the division rate (Salehi uses 1/day, the sensitivity analysis 0.141/day) affects the conclusion no less than β_hi.
4. **An alternative value:** for a more conservative upper limit, 3×10⁻³ could be used (the middle of the FISH range of CIN cell lines, ×60). The noise flag would then trigger more easily, but some genuinely high CIN would be misjudged as noise.

**Please confirm β_hi** (recommended 1×10⁻², or a value you specify); N1–N3 start after confirmation.

## References

1. Lengauer C, Kinzler KW, Vogelstein B. Nature 1997;386:623. doi:10.1038/386623a0
2. Thompson SL, Compton DA. J Cell Biol 2008;180:665. doi:10.1083/jcb.200712029
3. Bakhoum SF, Thompson SL, Manning AL, Compton DA. Nat Cell Biol 2009;11:27. doi:10.1038/ncb1809
4. Bakhoum SF et al. Clin Cancer Res 2011;17:7704. doi:10.1158/1078-0432.CCR-11-2049
5. Bolhaqueiro ACF et al. Nat Genet 2019;51:824. doi:10.1038/s41588-019-0399-6
6. Bollen Y et al. Nat Genet 2021;53:1187. doi:10.1038/s41588-021-00891-2
7. Lynch AR, Arp NL, Zhou AS, Weaver BA, Burkard ME. eLife 2022;11:e69799. doi:10.7554/eLife.69799
8. Lynch AR et al. PNAS 2024;121:e2309621121. doi:10.1073/pnas.2309621121
9. Janssen A, Kops GJPL, Medema RH. PNAS 2009;106:19108. doi:10.1073/pnas.0904343106
10. Silk AD et al. PNAS 2013;110:E4134. doi:10.1073/pnas.1317042110
11. Funk LC, Zasadil LM, Weaver BA. Dev Cell 2016;39:638. doi:10.1016/j.devcel.2016.10.023
12. Komarova NL, Wodarz D. PNAS 2004;101:7017. doi:10.1073/pnas.0401943101
13. Elizalde S, Laughney AM, Bakhoum SF. PLoS Comput Biol 2018;14:e1006447. doi:10.1371/journal.pcbi.1006447
14. Laughney AM et al. Cell Rep 2015;12:809. doi:10.1016/j.celrep.2015.06.065 (numbers not verified in the original)
15. Thompson SL, Compton DA. PNAS 2011;108:17974. doi:10.1073/pnas.1109720108
