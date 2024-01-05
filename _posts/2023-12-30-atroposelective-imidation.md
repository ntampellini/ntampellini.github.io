---
layout: post
title: "Atroposelective Imidation [ChemRXiv]"
tags: [interactive, computational_chemistry, publication]
---

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
    width: 50vmax;
    height: 30vmax;
    margin-left:auto;
    margin-right:auto;
  }

  /* Expand to the entire container */
  .iframe{
    width: 100%;
    height: 100%;
  }
</style>

<base href="/assets/atroposelective_imidation/" target="_blank">

<!-- CONTENT -->

Over the last ten years, we surely have become more used to seeing drug molecules featuring elements of axial chirality, with more and more of such drugs being developed as single stereoisomers. Still, if one axis of chirality is not something that turns heads anymore, a drug molecule featuring _two_ is sure to be a different story. Particularly if both these elements are precisely controlled during the process route and the drug is developed as a single, stable stereoisomer. This is the case for BMS-986142, an investigational BTK (Bruton's Tyrosine Kinase) inhibitor from Bristol Myers Squibb.[^BMSBTK]

<div class="img1cm">
  <!-- <img style='height: 40%; width: 40%; object-fit: contain' src="BMS.png"/> -->
  <img style='height: 20rem; width: 20rem; object-fit: contain' src="BMS.png"/>
  <br>
  <div align="center"><i><b>BMS-986142</b> - Investigational BTK inhibitor (leukemia treatment)</i>.</div>
</div>

The remarkable process route[^BMSBTK] starts by installing the <span class="ured">top right chiral center</span> via a perfectly enantioselective rhodium-catalyzed conjugate addition. Then, the <span class="uorange">biaryl chiral axis</span> is installed by another organometallic reaction with a chiral catalyst, namely a palladium-catalyzed Suzuki-Miyaura reaction with (_R_)-BINAP, achieving very high diastereoselectivity after extensive optimization. Interestingly, the <span class="ugreen">third chiral element</span> (second axis of chirality) is installed last, through a finely-orchestrated diastereoselective cyclization relying on the biaryl fluoride moiety. The authors justify this ordering of events with the relative stability of the two chiral axes, but we could not help noticing that atroposelective imidations to form quinazolinedione scaffolds were unknown at the time, except for a partially successful attempt from another BMS group at a Ni(0) catalyzed transformation that met with numerous difficulties.[^BMSNi0]

This is when we set out to develop an atroposelective, organocatalytic method to confer complete catalyst control over the formation of this novel scaffold. Inspired by the process route, we wanted to adopt a chiral superbase to affect the cyclization of a benzamide carbamate to a quinazolinedione. Our choice of a catalytic moiety eventually converged on tetramethylguanidines.

<div class="img1cm">
  <img style='height: 95%; width: 95%; object-fit: contain' src="strategy.png"/>
  <div align="center"><i>Atroposelective imidation strategy.</i></div>
</div>

<!-- As my first Ph.D. project in the lab of Prof. Miller, we developed an atroposelective imidation catalyst able to provide axially chiral quinazolinediones. The manuscript is currently in review and can be read on ChemRXiv[^chemrxiv]. Presented here is a summary of the highlights of the work along with some interactive renderings of the key geometries. Please refer to the full text for a comprehensive discussion. -->

A lot of work later, we were able to obtain a selective and specific catalyst, featuring novel and unusual key elements. Here is a sketch of the transformation and the lead catalyst.

<div class="img1cm">
  <img style='height: 95%; width: 95%; object-fit: contain' src="reaction.png"/>
</div>

Intrigued by the catalyst novelty, specificity, and the multi-step nature of the reaction, we decided to take a deep dive into the reaction mechanism with the help of DFT. First, we uncovered a conformational equilibrium of the catalyst between a _hairpin_ and a _folded_ form. While the former is the canonical β-turn, β-hairpin conformation, featuring a roughly planar arrangement of the four amino acid residues, the latter _folded_ conformation features an additional internal hydrogen bond, conferring the structure a helical, three-dimensional character. More importantly, the _hairpin_ foldamer features a hydrogen bond between the guanidine and the difluoroacetyl group on proline, while the _folded_ state does not. We believe that this effect explains why previous generations of the catalyst were suddenly completely inactive: they were heavily favoring the _hairpin_ conformer, where the active guanidine residue was masked and not available for catalysis. Below you can find an animation connecting the two aforementioned conformers.

> Click and drag to rotate, Ctrl+click/Mouse3 to move, scroll/right-click to zoom.

### Catalyst conformational dynamism

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="py3Dmol/folding.html"
    frameBorder="0"
  ></iframe>
</div>

<div align="center"><i>Lead catalyst converting between its hairpin and folded conformations.</i></div>

Moving on to the transformation, we devised a large reaction space to be characterized computationally, and we ended up modeling the two reaction steps (addition and elimination) through ten different activation modes across four substrate stereoisomers and two catalyst foldamers, for a total of ~400 transition states. After achieving a solid alignment with the experimental results, we started looking for the selectivity-defining factors. An initial, striking observation we made was that the four productive transition states leading to the two enantiomers of the product **all featured topologically different transition states**. Each enantiomer of the product is formed via a different diastereoisomer and features a unique activation mode for each of the two reaction steps. Moreover, each enantiomer of the product also has a different rate-determining step!

<div class="img1cm">
  <img style='height: 100%; width: 100%; object-fit: contain' src="PES.png"/>

  <div align="center"><i>Productive potential energy surface of the reaction. Energetic data at the M06-2X/def2-QZVP/CPCM(PhCF<sub>3</sub>)//R<sup>2</sup>SCAN-3c/CPCM(PhCF<sub>3</sub>) level.</i></div>
</div>

How can we tease apart individual effects in such a complex reaction landscape? By closely inspecting transition state structures, we noticed the importance of the conformation of the unusual difluoroacetyl moiety present in the catalyst. The close contacts between the strongly polarized C(F<sub>2</sub>)-H bond and oxygen atoms hinted at the presence of CH-O non-classical hydrogen bonds (2.7-3.1 Å). This group showed the remarkable ability to selectively engage as a _bidentate_ hydrogen bond donor **only** in the transition states leading to the _major_ enantiomer (S<sub>a</sub>), while only acting as _monodentate_ in the disfavored pathway to the _minor_ (R<sub>a</sub>). Here you can find some interactive models of the transition states, where you can observe this behavior.

### TS1 S<sub>a</sub> (major, bidentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="py3Dmol/TS1Sa.html"
    frameBorder="0"
  ></iframe>
</div>

### TS1 R<sub>a</sub> (minor, monodentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="py3Dmol/TS1Ra.html"
    frameBorder="0"
  ></iframe>
</div>

In the second transition state (TS2, _elimination_) the situation is the same, and the _major_ enantiomer (S<sub>a</sub>) is formed via an activation mode that features bidentate coordination of the difluoroacetamide moiety. The _minor_ enantiomer once again only features monodentate coordination in the transition state (TS2 R<sub>a</sub>).

### TS2 S<sub>a</sub> (major, bidentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="py3Dmol/TS2Sa.html"
    frameBorder="0"
  ></iframe>
</div>

### TS2 R<sub>a</sub> (minor, monodentate difluoroacetamide)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="py3Dmol/TS2Ra.html"
    frameBorder="0"
  ></iframe>
</div>

All the structures can be found in the [GitHub repository](https://github.com/ntampellini/Atroposelective_Imidation) of this work.

We are currently leveraging the conformational insights and the catalyst design we developed in other asymmetric reactions.

---
{: data-content="footnotes"}

[^BMSBTK]: [Discovery](https://pubs.acs.org/doi/full/10.1021/acs.jmedchem.6b01088), [synthesis](https://pubs.acs.org/doi/full/10.1021/acs.orglett.8b01218) and [process optimization](https://pubs.acs.org/doi/10.1021/acs.oprd.8b00246) of BMS-986142.

[^BMSNi0]: [Ni(0)-catalyzed synthesis of quinazolinediones](https://pubs.acs.org/doi/full/10.1021/acs.orglett.7b00052).

[^chemrxiv]: You can find the preprint on [ChemRXiv](https://chemrxiv.org/engage/chemrxiv/article-details/658102309138d231610176ba).

