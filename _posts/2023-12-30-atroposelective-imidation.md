---
layout: post
title: "Atroposelective Imidation [ChemRXiv]"
tags: [interactive, computational_chemistry, publication]
---

As my first Ph.D. project in the lab of Prof. Miller, we developed an atroposelective imidation catalyst able to provide axially chiral quinazolinediones. The manuscript is currently in review and can be read on ChemRXiv[^chemrxiv]. Presented here is a summary of the highlights of the work along with some interactive renderings of the key geometries. Please refer to the full text for a comprehensive discussion.

Following a series of structure-activity relationships we observed experimentally, we obtained a selective and specific catalyst featuring novel key elements. Here is a sketch of the transformation.

<div>
  <img style='height: 100%; width: 100%; object-fit: contain' src="/assets/atroposelective_imidation/reaction.png"/>
</div>

Intrigued by the catalyst novelty, specificity, and the multi-step nature of the reaction, we decided to take a deep dive into the reaction mechanism with the help of DFT. First, we uncovered a conformational equilibrium of the catalyst between a _hairpin_ and a _folded_ form. While the former is the canonical β-turn, β-hairpin conformation, featuring a roughly planar arrangement of the four amino acid residues, the latter _folded_ conformation features an additional internal hydrogen bond, conferring the structure a helical, three-dimensional character. More importantly, the _hairpin_ foldamer features a hydrogen bond between the guanidine and the difluoroacetyl group on proline, while the _folded_ state does not. We believe that this effect explains why previous generations of the catalyst were suddenly completely inactive: they were heavily favoring the _hairpin_ conformer, where the active guanidine residue was masked and not available for catalysis. Below you can find an animation connecting the two aforementioned conformers.

> Click and drag to rotate, Ctrl+click to move, scroll to zoom.

### Catalyst conformational dynamism
<iframe src="/assets/atroposelective_imidation/py3Dmol/folding.html" width="675" height="425" frameBorder="0"></iframe>
<!-- <iframe src="/assets/atroposelective_imidation/py3Dmol/catalysts_static.html" width="675" height="425" frameBorder="0"></iframe> -->

<div align="center"><i>Lead catalyst in the hairpin (left) and folded (right) conformations.</i></div>

Moving on to the transformation, we devised a large reaction space to be characterized computationally, and we ended up modeling the two reaction steps (addition and elimination) through ten different activation modes across four substrate stereoisomers and two catalyst foldamers, for a total of ~400 transition states. After achieving a solid alignment with the experimental results, we started looking for the selectivity-defining factors. An initial, striking observation we made was that the four productive transition states leading to the two enantiomers of the product **all featured topologically different transition states**. Each enantiomer of the product is formed via a different diastereoisomer and features a unique activation mode for each of the two reaction steps. Moreover, each enantiomer of the product also has a different rate-determining step!

<div>
  <img style='height: 100%; width: 100%; object-fit: contain' src="/assets/atroposelective_imidation/PES.png"/>
</div>

<div align="center"><i>Productive potential energy surface of the reaction. Energetic data at the M06-2X/def2-QZVP/CPCM(PhCF<sub>3</sub>)//R<sup>2</sup>SCAN-3c/CPCM(PhCF<sub>3</sub>) level.</i></div>

How can we tease apart individual effects in such a complex reaction landscape? By closely inspecting transition state structures, we noticed the importance of the conformation of the unusual difluoroacetyl moiety present in the catalyst. The close contacts between the strongly polarized C(F<sub>2</sub>)-H bond and oxygen atoms hinted at the presence of CH-O non-classical hydrogen bonds (2.7-3.1 Å). This group showed the remarkable ability to selectively engage as a _bidentate_ hydrogen bond donor **only** in the transition states leading to the _major_ enantiomer (S<sub>a</sub>), while only acting as _monodentate_ in the disfavored pathway to the _minor_ (R<sub>a</sub>). Here you can find some interactive models of the transition states, where you can observe this behavior.

### TS1 S<sub>a</sub> (major, bidentate difluoroacetamide)
<iframe src="/assets/atroposelective_imidation/py3Dmol/TS1Sa.html" width="675" height="425" frameBorder="0"></iframe>

### TS1 R<sub>a</sub> (minor, monodentate difluoroacetamide)
<iframe src="/assets/atroposelective_imidation/py3Dmol/TS1Ra.html" width="675" height="425" frameBorder="0"></iframe>

In the second transition state (TS2, _elimination_) the situation is the same, and the _major_ enantiomer (S<sub>a</sub>) is formed via an activation mode that features bidentate coordination of the difluoroacetamide moiety. The _minor_ enantiomer once again only features monodentate coordination in the transition state (TS2 R<sub>a</sub>).

### TS2 S<sub>a</sub> (major, bidentate difluoroacetamide)
<iframe src="/assets/atroposelective_imidation/py3Dmol/TS2Sa.html" width="675" height="425" frameBorder="0"></iframe>

### TS2 R<sub>a</sub> (minor, monodentate difluoroacetamide)
<iframe src="/assets/atroposelective_imidation/py3Dmol/TS2Ra.html" width="675" height="425" frameBorder="0"></iframe>

All the structures can be found in the [GitHub repository](https://github.com/ntampellini/Atroposelective_Imidation) of this work.

We are currently leveraging the conformational insights and the catalyst design we developed in other asymmetric reactions.

---
{: data-content="footnotes"}

[^chemrxiv]: You can find the preprint on [ChemRXiv](https://chemrxiv.org/engage/chemrxiv/article-details/658102309138d231610176ba).

