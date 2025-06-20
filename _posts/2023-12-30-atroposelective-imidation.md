---
layout: post
title: "Atroposelective Imidation of Quinazolinediones"
tags: [interactive, computational_chemistry, publication]
---

<!-- Twitter -->
<meta property="twitter:card" content="summary_large_image" />
<meta property="twitter:url" content="https://ntampellini.github.io/atroposelective-imidation/" />
<meta property="twitter:title" content="Catalyst Conformational Activation: Atroposelective Imidation of Quinazolinediones" />
<meta property="twitter:description" content="Test description" />
<meta property="twitter:image" content="/assets/atroposelective_imidation/TS1Sa_preview.png" />

<!-- Meta Tags Generated with https://metatags.io -->

<style>
  .img1cm {
     margin-top: 1cm;
     margin-bottom: 1cm;
  }

  .ured {
    font-weight: bold;
    text-decoration: underline;
    text-decoration-color: rgba(232, 111, 136, 0.75);
    text-decoration-thickness: 0.25em;
    text-decoration-skip-ink: none;
  }

  .uorange {
    font-weight: bold;
    text-decoration: underline;
    text-decoration-color: rgba(255, 195, 135, 0.75);
    text-decoration-thickness: 0.25em;
    text-decoration-skip-ink: none;
  }

  .ugreen {
    font-weight: bold;
    text-decoration: underline;
    text-decoration-color: rgba(100, 200, 200, 0.75);
    text-decoration-thickness: 0.25em;
    text-decoration-skip-ink: none;
  }

  .iframe-parent{
    max-width: min(90vw, 1000px);
    max-height: min(80vh, 700px);
    width: 95vw;
    height: 75vw;
    margin-left:auto;
    margin-right:auto;
  }

  /* Expand to the entire container */
  .iframe{
    width: 100%;
    height: 100%;
  }
</style>

<!-- CONTENT -->

> This post summarizes and adds interactive visualizations on the development of a novel superbasic peptidyl guanidine catalyst able to affect the atroposelective imidation of quinazolinediones. The full paper, previously on ChemRXiv[^chemrxiv] is now published in Chemistry, a European Journal.[^ChemEurJ]

> Click and drag to rotate, Ctrl+click/Mouse3 to move, scroll/right-click to zoom.

Over the last ten years, we surely have become more used to seeing drug molecules featuring elements of axial chirality. Still, if a chiral axis is not something that turns heads anymore, a drug molecule featuring _two_ is sure to be a different story. Particularly if both these elements are precisely controlled during the process route and the drug is developed as a single, stable stereoisomer. This is the case for BMS-986142, an investigational BTK (Bruton's Tyrosine Kinase) inhibitor from Bristol Myers Squibb.[^BMSBTK]

<div class="img1cm">
  <!-- <img style='height: 40%; width: 40%; object-fit: contain' src="BMS.png"/> -->
  <img style='max-width: 95vw; width: 300px; object-fit: contain' src="/assets/atroposelective_imidation/BMS.png"/>
  <br>
  <div align="center"><i><b>BMS-986142</b> - Investigational BTK inhibitor (leukemia treatment)</i>.</div>
</div>

The remarkable process route[^BMSBTK] starts by installing the <span class="ured">top right chiral center</span> via a perfectly enantioselective rhodium-catalyzed conjugate addition. Then, the <span class="uorange">biaryl chiral axis</span> is installed by another organometallic reaction with a chiral catalyst, namely a palladium-catalyzed Suzuki-Miyaura reaction with (_R_)-BINAP, achieving very high diastereoselectivity. Interestingly, the <span class="ugreen">third chiral element</span> (_i.e._ second axis of chirality) is installed last, through a finely-orchestrated diastereoselective cyclization relying on the biaryl fluoride moiety. The authors explain this ordering with the relative stability of the two chiral axes, but we could not help noticing that atroposelective imidations to form quinazolinedione scaffolds were almost unknown at the time. Indeed, the only exception is a partially successful attempt from another BMS group, adopting a Ni(0)-catalyzed isocyanate addition that met with various difficulties, particularly in the tentative enantioselective extension.[^BMSNi0]

This perspective set us out to seek a different solution for this problem, given the now increasing presence of chiral axes in medicinally relevant, heterocyclic scaffolds. Our vision was to develop an atroposelective, organocatalytic method to confer complete catalyst control over the formation of the C-N chiral axis. Inspired by BMS-986142's process route, we inquired on adopting a chiral superbase to affect the cyclization of a benzamide carbamate to a quinazolinedione. Our choice of a catalytic moiety eventually converged on tetramethylguanidines.

<div class="img1cm">
  <img style='max-width: 95vw; width: 900px; object-fit: contain' src="/assets/atroposelective_imidation/strategy.png"/>
  <div align="center"><i>Atroposelective imidation strategy.</i></div>
</div>

Following extensive catalyst optimization, we were able to obtain a selective and specific catalyst, featuring novel and unusual key elements. Here is a sketch of the transformation and the lead catalyst.

<div class="img1cm">
  <img style='max-width: 95vw; width: 900px; object-fit: contain' src="/assets/atroposelective_imidation/reaction.png"/>
</div>

Intrigued by the catalyst novelty, specificity, and the multi-step nature of the reaction, we decided to take a deep dive into the reaction mechanism with the help of DFT. First, we uncovered a conformational equilibrium of the catalyst between a _hairpin_ and a _folded_ form. While the former is the canonical β-turn, β-hairpin conformation, featuring a roughly planar arrangement of the four amino acid residues, the latter _folded_ conformation features an additional internal hydrogen bond, conferring the structure a helical, three-dimensional character. More importantly, the _hairpin_ foldamer features a hydrogen bond between the guanidine and the difluoroacetyl group on proline, while the _folded_ state does not. We believe that this effect explains why previous generations of the catalyst were suddenly completely inactive: they were heavily favoring the _hairpin_ conformer, where the active guanidine residue was masked and not available for catalysis. Below you can find an animation connecting the two aforementioned conformers.

> Click and drag to rotate, Ctrl+click/Mouse3 to move, scroll/right-click to zoom.

### Catalyst conformational dynamism

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/atroposelective_imidation/py3Dmol/folding.html"
    frameBorder="0"
  ></iframe>
</div>

<div align="center"><i>Lead catalyst converting between its hairpin and folded conformations.</i></div>

Moving on to the transformation, we devised a large reaction space to be characterized computationally, and we ended up modeling the two reaction steps (addition and elimination) through ten different activation modes across four substrate stereoisomers and two catalyst foldamers, for a total of ~400 transition states. After achieving a solid alignment with the experimental results, we started looking for the selectivity-defining factors. An initial, striking observation we made was that the four productive transition states leading to the two enantiomers of the product **all featured topologically different transition states**. Each enantiomer of the product is formed via a different diastereoisomer and features a unique activation mode for each of the two reaction steps. Moreover, each enantiomer of the product also has a different rate-determining step!

<div class="img1cm">
  <img style='max-width: 95vw; width: 900px; object-fit: contain' src="/assets/atroposelective_imidation/PES.png"/>

  <div align="center"><i>Productive potential energy surface of the reaction. Energetic data at the M06-2X/def2-QZVP/CPCM(PhCF<sub>3</sub>)//R<sup>2</sup>SCAN-3c/CPCM(PhCF<sub>3</sub>) level.</i></div>
</div>

How can we tease apart individual effects in such a complex reaction landscape? By closely inspecting transition state structures, we noticed the importance of the conformation of the unusual difluoroacetyl moiety present in the catalyst. The close contacts between the strongly polarized C(F<sub>2</sub>)-H bond and oxygen atoms hinted at the presence of CH-O non-classical hydrogen bonds (2.7-3.1 Å). This group showed the remarkable ability to selectively engage as a _bidentate_ hydrogen bond donor **only** in the transition states leading to the _major_ enantiomer (S<sub>a</sub>), while only acting as _monodentate_ in the disfavored pathway to the _minor_ (R<sub>a</sub>). Here you can find some interactive models of the transition states, where you can observe this behavior.

### TS1 S<sub>a</sub> (major, bidentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/atroposelective_imidation/py3Dmol/TS1Sa.html"
    frameBorder="0"
  ></iframe>
</div>

### TS1 R<sub>a</sub> (minor, monodentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/atroposelective_imidation/py3Dmol/TS1Ra.html"
    frameBorder="0"
  ></iframe>
</div>

In the second transition state (TS2, _elimination_) the situation is the same, and the _major_ enantiomer (S<sub>a</sub>) is formed via an activation mode that features bidentate coordination of the difluoroacetamide moiety. The _minor_ enantiomer once again only features monodentate coordination in the transition state (TS2 R<sub>a</sub>).

### TS2 S<sub>a</sub> (major, bidentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/atroposelective_imidation/py3Dmol/TS2Sa.html"
    frameBorder="0"
  ></iframe>
</div>

### TS2 R<sub>a</sub> (minor, monodentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/atroposelective_imidation/py3Dmol/TS2Ra.html"
    frameBorder="0"
  ></iframe>
</div>

### Conclusions

After learing so much on this transformation and its relationship with the catalyst, we are eager to initiate the knowledge transfer on other systems. We are currently leveraging these conformational insights to guide the development of new chiral catalysts, tailored to enable new, selective reactions. By understanding the key elements dictating selectivity, we were able to dissect complex relationships that would be very hard to parse out from the experimental data alone. It is my hope that understanding the conformational dynamism of this new catalyst will pave the way to even deeper insights on these peptide catalysts, and will more broadly help to speed up the ever-important process of catalyst optimization, that often feels like _looking for needle(s) in a peptide haystack._

All the structures in this post can be found in the [GitHub repository](https://github.com/ntampellini/Atroposelective_Imidation) of this work.

---
{: data-content="footnotes"}

[^BMSBTK]: [Discovery](https://pubs.acs.org/doi/full/10.1021/acs.jmedchem.6b01088), [synthesis](https://pubs.acs.org/doi/full/10.1021/acs.orglett.8b01218) and [process optimization](https://pubs.acs.org/doi/10.1021/acs.oprd.8b00246) of BMS-986142.

[^BMSNi0]: [Ni(0)-catalyzed synthesis of quinazolinediones](https://pubs.acs.org/doi/full/10.1021/acs.orglett.7b00052).

[^chemrxiv]: The preprint of this work can be found on [ChemRXiv](https://chemrxiv.org/engage/chemrxiv/article-details/658102309138d231610176ba).

[^ChemEurJ]: The full paper can be found online on [Chemistry, a European Journal](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.202401109).