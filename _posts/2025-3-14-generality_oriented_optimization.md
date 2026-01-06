---
layout: post
title: "Computational Analysis of a Generality-Oriented Optimization Campaign"
tags: [interactive, computational_chemistry, publication]
---

<!-- Twitter -->
<!-- <meta property="twitter:card" content="summary_large_image" />
<meta property="twitter:url" content="https://ntampellini.github.io/enantioselective_s_oxidation/" />
<meta property="twitter:title" content="Enantioselective S-Oxidation with Aspartic acid-derived peptides" />
<meta property="twitter:description" content="Test description" /> -->
<!-- <meta property="twitter:image" content="/assets/atroposelective_imidation/TS1Sa_preview.png" />

<!-- Meta Tags Generated with https://metatags.io -->

<style>
  .img1cm {
     margin-top: 1cm;
     margin-bottom: 1cm;
  }

  .chemdraw {
    width: 85vw;
    max-width: 800px;
    max-height: 90vh;
    object-fit: contain;
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

  hr{
    border: 0;
    border-bottom: 2px solid rgba(75, 75, 75, 255);
    height: 3px;
  }

.tab-buttons {
  display: flex;
}

.tab-button {
  background-color: var(--main-back-color);
  border: none;
  border-radius: 5px;
  border-bottom-left-radius: 0px;
  border-bottom-right-radius: 0px;
  padding: 10px 20px;
  cursor: pointer;
  height: 40px;
}

.tab-button:hover {
  /* background-color: var(--target-highlight-color); */
}

.tab {
  display: none;
  border-radius: 5px;
  border-top-left-radius: 0px;
  width: 95%;
  padding: 20px;
  background-color: var(--quote-background);
}

.active-tab-button {
  background-color: var(--quote-background);
  font-weight: bold;
  /* transform: scaleY(1.1); */
  /* transform-origin: bottom; */
}

/* Claude mobile-friendly CSS layout only on mobile */

/* Add smooth horizontal scrolling on mobile */
.tab-buttons {
  display: flex;
  overflow-x: auto;
  overflow-y: hidden;
  -webkit-overflow-scrolling: touch; /* smooth scrolling on iOS */
  scrollbar-width: thin; /* Firefox */
  gap: 2px; /* small space between buttons */
}

/* Hide scrollbar on Chrome/Safari but keep functionality */
.tab-buttons::-webkit-scrollbar {
  height: 4px;
}

.tab-buttons::-webkit-scrollbar-track {
  background: transparent;
}

.tab-buttons::-webkit-scrollbar-thumb {
  background: rgba(128, 128, 128, 0.3);
  border-radius: 2px;
}

.tab-button {
  background-color: var(--main-back-color);
  border: none;
  border-radius: 5px;
  border-bottom-left-radius: 0px;
  border-bottom-right-radius: 0px;
  padding: 10px 20px;
  cursor: pointer;
  height: 40px;
  flex-shrink: 0;
  white-space: nowrap;
}

/* Placing this AFTER .tab-button to ensure it takes precedence */
.tab-button.active-tab-button {
  background-color: var(--quote-background);
  font-weight: bold;
}

/* reduce padding on very small screens */
@media (max-width: 480px) {
  .tab-button {
    padding: 10px 12px;
    font-size: 0.9em;
  }
}

</style>

<!-- CONTENT -->

> This post reuses figures and tables from an ACS Catalysis publication[^ref]. Reproduction follows the ACS guidelines.[^ACS_guidelines]

Here you can find interactive visualizations for our computational analysis on the "Experimental Lineage and Computational Analysis of a General Aminoxyl-based Oxidation Catalyst: Generality from Substrate-Specific Interactions"[^ref] by Rozema, S. D.<sup>‡</sup>; Tampellini, N.<sup>‡</sup>; Rein, J.; Sigman, M. S.; Lin, S. and Miller, S. J.*.

## Single point changes: from P1 to P7

<div id="panel0" class="tab-buttons">
    <button class="tab-button active-tab-button" onclick="toggleTab('panel0', 0)">Introduction</button>
    <button class="tab-button" onclick="toggleTab('panel0', 1)">P1 → P3</button>
    <button class="tab-button" onclick="toggleTab('panel0', 2)">P2 → P4</button>
    <button class="tab-button" onclick="toggleTab('panel0', 3)">P3 → P4</button>
    <button class="tab-button" onclick="toggleTab('panel0', 4)">P4 → P5</button>
    <button class="tab-button" onclick="toggleTab('panel0', 5)">P6 → P7</button>
</div>

<div class="tab" id="tab1">
  <p style="margin-top: 0px"> The ground state interaction mode of pairs of peptide catalysts were analyzed for a single substrate. A substrate was chosen for each peptide pair such that the generational advancement of the catalyst would induce the greatest improvement in selectivity for that substrate (and other considerations, see main text). Then, peptide-substrate covalent adducts were modeled as four possible stereoisomers, as alcohol binding to the oxoammonium moiety generates two new chiral elements: the first arising from the choice of which of the two meso alcohol groups is bound to the oxoammonium moiety, and the second from which oxoammonium face was used for coordination. The >6k structures modeling at DFT level and their conformational analysis helped uncover the atomistic rationales underpinning the success of the catalyst optimization campaign, encompassing an overall trend of conformational complexity reduction and convergence towards to a single activation mode. Yet, not all selective substrates interact in the same way (<i>i.e.</i> with the same binding topology)!
  </p>
  <img style='width: 90vw; max-width: 750px; object-fit: contain' src="/assets/generality_oriented_optimization/table.jpg"/>
</div>

<div class="tab" id="tab2">
  <h3 style="margin-top: 0px;"> P1 → P3 (NMe<sub>2</sub> → NHMe) </h3>
  <div>
    <img class="chemdraw" src="/assets/generality_oriented_optimization/p1p3.jpg"/>
    <!-- <div align="center"><i><b>P1 to P3</b> - Description</i>.</div> -->
  </div>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p1p3.html"
      frameBorder="0"
    ></iframe>
  </div>
  <div align="center"><i>Mouse to rotate, Ctrl+click to move, scroll/right-click to zoom.</i></div>
</div>

<div class="tab" id="tab3">
  <h3 style="margin-top: 0px;">P2 → P4 (NMe<sub>2</sub> → NHMe)</h3>
  <div>
    <img class="chemdraw" src="/assets/generality_oriented_optimization/p2p4.jpg"/>
    <!-- <div align="center"><i><b>P2 to P4</b> - Description</i>.</div> -->
  </div>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p2p4.html"
      frameBorder="0"
    ></iframe>
  </div>
  <div align="center"><i>Mouse to rotate, Ctrl+click to move, scroll/right-click to zoom.</i></div>
</div>

<div class="tab" id="tab4">
  <h3 style="margin-top: 0px;">P3 → P4 (<sup>D</sup>Phe → <sup>L</sup>Phe)</h3>
  <div>
    <img class="chemdraw" src="/assets/generality_oriented_optimization/p3p4.jpg"/>
    <!-- <div align="center"><i><b>P3 to P4</b> - Description</i>.</div> -->
  </div>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p3p4.html"
      frameBorder="0"
    ></iframe>
  </div>
  <div align="center"><i>Mouse to rotate, Ctrl+click to move, scroll/right-click to zoom.</i></div>
</div>

<div class="tab" id="tab5">
  <h3 style="margin-top: 0px;">P4 → P5 (Pro → Pip)</h3>
  <div>
    <img class="chemdraw" src="/assets/generality_oriented_optimization/p4p5.jpg"/>
    <!-- <div align="center"><i><b>P1 to P3</b> - Description</i>.</div> -->
  </div>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p4p5.html"
      frameBorder="0"
    ></iframe>
  </div>
  <div align="center"><i>Mouse to rotate, Ctrl+click to move, scroll/right-click to zoom.</i></div>
</div>

<div class="tab" id="tab6">
  <h3 style="margin-top: 0px;">P5 → P7 (Phe → Bip)</h3>
  <div>
    <img class="chemdraw" src="/assets/generality_oriented_optimization/p5p7.jpg"/>
    <!-- <div align="center"><i><b>P5 to P7</b> - Description</i>.</div> -->
  </div>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p5p7.html"
      frameBorder="0"
    ></iframe>
  </div>
  <div align="center"><i>Mouse to rotate, Ctrl+click to move, scroll/right-click to zoom.</i></div>
</div>

## Generality from substrate-specific interactions

<div id="panel1" class="tab-buttons">
    <button class="tab-button active-tab-button" onclick="toggleTab('panel1', 0)">Why is <b>P7</b> so general?</button>
    <button class="tab-button" onclick="toggleTab('panel1', 1)">S4 (96% ee)</button>
    <button class="tab-button" onclick="toggleTab('panel1', 2)">S10 (90% ee)</button>
    <button class="tab-button" onclick="toggleTab('panel1', 3)">S15 (26% ee)</button>
</div>

<div class="tab" id="tab1">
  <p style="margin-top: 0px">After observing the convergence of many substrates to a unified interaction topology across the different catalyst generations, we wanted to assess if multiple substrates that are effectively desymmetrized by the optimized catalyst <b>P7</b> also favor the same interaction network, in addition to the already modeled <b>S12</b>. Surprisingly, the situation proved to be more nuanced, with different substrates showing different interaction modes. The generality of high selectivity for many substrates was found to be consistent with alternative, and indeed substrate-specific interactions, suggesting that functional generality need not be a result of strict mechanistic homology.</p>
  <div>
    <img class="chemdraw" style="margin-bottom: 10px" src="/assets/generality_oriented_optimization/modes.jpg"/>
    <div align="center"><i>Substrate generality was achieved also through (or despite) substrate-specific interactions, two mechanistic features often assumed to be incompatible.</i></div>
  </div>
</div>

<div class="tab" id="tab2">
  <h3 style="margin-top: 0px;">P7 / S4 (96% ee)</h3>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p7s4.html"
      frameBorder="0"
    ></iframe>
  </div>
</div>

<div class="tab" id="tab3">
  <h3 style="margin-top: 0px;">P7 / S10 (90% ee)</h3>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p7s10.html"
      frameBorder="0"
    ></iframe>
  </div>
</div>

<div class="tab" id="tab4">
  <h3 style="margin-top: 0px;">P7 / S15 (26% ee)</h3>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/generality_oriented_optimization/p7s15.html"
      frameBorder="0"
    ></iframe>
  </div>
</div>

### Extension to regioselective primary alcohol oxidation

The Lin and Miller groups also reported the use of an oxo-analog of <b>P7</b> (<b>OAzcP7</b>) to enable the regiocontrolled oxidation of unsymmetrical diols[^regio]. In this reaction, the highly preorganized peptidyl template enabled the selective oxidation of the more hindered alcohol group over the more accessible one, overriding intrinsic steric biases. Ground state modeling of this system demonstrated its operation through conceptually analogous preferential formation of one catalyst-substrate adduct over its regioisomer.

<div>
  <img class="chemdraw" src="/assets/generality_oriented_optimization/p7d2.jpg"/>
  <!-- <div align="center"><i>The same optimized catalyst <b>P7</b> can also successfuly tackle the regioselective oxidation of <b>D2</b></i></div> -->
</div>

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/generality_oriented_optimization/p7d2.html"
    frameBorder="0"
  ></iframe>
</div>

<script src="/assets/tabs.js"></script>

---
{: data-content="footnotes"}

[^ref]: Read the full paper on [ACS Catalysis](https://pubs.acs.org/doi/10.1021/acscatal.5c05893)
[^ACS_guidelines]: See [ACS guidelines for scholarly posting & sharing policies](https://pubs.acs.org/page/copyright/journals/posting_policies.html)
[^regio]: Catalyst-Controlled Regiodivergent Oxidation of Unsymmetrical Diols, read at [JACS](https://pubs.acs.org/doi/10.1021/jacs.5c00330)