---
layout: post
title: "Catalytic Asymmetric Amination of Sulfenamides by DFT-guided Catalyst Design"
tags: [interactive, computational_chemistry, publication]
---

<style>
  .img1cm {
     margin-top: 1cm;
     margin-bottom: 1cm;
  }

  .chemdraw {
    width: 85vw;
    max-width: 1000px;
    max-height: 90vh;
    object-fit: contain;
    margin-bottom: 0.5cm;
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
    overflow: hidden;
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

  /* responsive video */
  .responsive-video {
    height: 0;
    padding-top: 1px;
    position: relative;
    margin-bottom: -3px;
    padding-bottom: 56.25%;
  }
  .responsive-video iframe {
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      position: absolute;
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

> This post reuses figures and tables from an ACS publication[^ref]. Reproduction follows the ACS guidelines.[^ACS_guidelines]

<!-- > Click and drag to rotate, Ctrl+click to move, scroll/right-click to zoom. -->

This blog post contains interactive visualization for the "Peptide-Catalyzed Asymmetric Amination of Sulfenamides Enabled by DFT-Guided Catalyst Optimization"[^ref] by Tampellini, N.; Choi, E. S. and Miller, S. J.* (<i>J. Am. Chem. Soc.</i> <b>2025</b>)[^ref].

<div>
  <img class="chemdraw" src="/assets/electrophilic_s_amination/reaction.jpg"/>
  <div align="center" style="margin-bottom: 0.5cm"><i>The electrophilic amination reaction of this work.</i></div>
</div>

<div id="panel0" class="tab-buttons">
    <button class="tab-button active-tab-button" onclick="toggleTab('panel0', 0)">Reaction concept</button>
    <button class="tab-button" onclick="toggleTab('panel0', 1)">Catalyst screen</button>
    <button class="tab-button" onclick="toggleTab('panel0', 2)">DFT: increase ΔΔG<sup>‡</sup></button>
    <button class="tab-button" onclick="toggleTab('panel0', 3)">DFT: reduce ΔG°<sub>turnover</sub></button>
</div>

<div class="tab" id="tab1">
  <p style="margin-top: 0px">Recent times have seen a new interest in chiral-at-sulfur functional groups, but asymmetric methods to prepare them remain scarce. After a report from Saito and coworkers on the electrophilic amination of sulfenamides with a chiral reagent, we realized the lack of a catalytic, asymmetric approach to sulfinamidines. Instead of acidic activation, we reasoned that the reaction could proceed under base catalysis as well. The energetically raised, diffused HOMO of the sulfenamidate anion would therefore react with an electrophilic aminating reagent. Gladly, reactivity was straightforward to establish - but enantioselectivity proved more challenging.
  </p>
  <img class="chemdraw" src="/assets/electrophilic_s_amination/figure1.jpg"/>
</div>

<div class="tab" id="tab2">
  <p style="margin-top: 0px">The screen of established chiral bases to direct the enantioselective amination was initially not very promising: only a few, mild hits emerged from a sensibly sized library. While peptidyl guanidines appeared as the most promising scaffold, our diverse library showed low variability in terms of enantioselectivity, with all catalysts performing similarly (~0-30% ee). We clearly didn't have the right catalysts for this transformation, and perhaps the community also did not. So how do we get to one that works? We challenged ourselves going for a DFT-first optimization approach, seeking to understand what lies beneath the mild (and unique) success of peptidyl guanidines as well as what doesn't work very well. From there, we can have a better idea of the catalyst chemical space we would need to explore.
  </p>
  <img class="chemdraw" src="/assets/electrophilic_s_amination/figure2.jpg"/>
</div>

<div class="tab" id="tab3">
  <p style="margin-top: 0px">Hit catalyst <b>4</b>, as expected, featured many energetically accessible activation modes - we targeted one as the target of our computational optimization. This mode positioned groups in a way that we believed could accomodate numerous substrates, while directing the amination reaction. Hence, we tested the strategic addition of steric bulk to a pendant amide position, with the idea of limiting the conformational space of the catalyst-substrate complex to that single mode. This would increase the ΔΔG<sup>‡</sup> and ideally increase enantioselectivity. After just a few iterations, DFT calculations suggested that we had done it. Single activation mode, ΔΔG<sup>‡</sup> > 4 kcal/mol! But wait, why is the selectivity still only about 4:1 er, then?
  </p>
  <img class="chemdraw" src="/assets/electrophilic_s_amination/figure3.jpg"/>
</div>

<div class="tab" id="tab4">
  <p style="margin-top: 0px">It turns out that we had done <i>too much</i> of a good job at designing the catalyst. It proved, in fact, to be such a good binder of the transition state that it was also binding the stereoelectronically  similar products (post-reaction complex) too tightly. This is a problem if you are developing a catalytic reaction, where you need to displace the products from the catalyst to start a new cycle. DFT thermochemistry, in fact, informed us that catalyst turnover (exchanging the reaction products for the reagents) was endergonic for the substrates we were looking at. This meant that after some initial progress, the reaction products would inhibit the catalyst and other non-enantioselective processes would mediate the unselective formation of our chiral sulfinamidine products. The key to escape this issue was to find substrates that had an easier time dissociating from the catalyst. These turned out to be meta-substituted benzoyl sulfinamidines. While this seems a little specific, this moiety is easy to install and remove later in derivatizations - practically a protecting group!
  </p>
  <img class="chemdraw" src="/assets/electrophilic_s_amination/figure4.jpg"/>
</div>

<p style="margin-bottom: 0.5cm"></p>

<div id="panel1" class="tab-buttons">
    <button class="tab-button active-tab-button" onclick="toggleTab('panel1', 0)">Major product TS</button>
    <button class="tab-button" onclick="toggleTab('panel1', 1)">Minor product TS</button>
    <button class="tab-button" onclick="toggleTab('panel1', 2)">Animation</button>
    <!-- <button class="tab-button" onclick="toggleTab('panel1', 3)">DFT: reduce ΔG°<sub>turnover</sub></button> -->
</div>

<div class="tab" id="tab1">
  <h3 style="margin-top: 0px;">Major product TS (pro-<i>S</i>)</h3>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/electrophilic_s_amination/major.html"
    frameBorder="0"
    ></iframe>
  </div>
<div align="center"><i>Substrate and aminating reagent highlighted in orange.</i></div>
</div>

<div class="tab" id="tab2">
  <h3 style="margin-top: 0px;">Minor product TS (pro-<i>R</i>)</h3>
  <div class="iframe-parent">
    <iframe
      class="iframe"
      src="/assets/electrophilic_s_amination/minor.html"
    frameBorder="0"
    ></iframe>
  </div>
<div align="center"><i>Substrate and aminating reagent highlighted in orange.</i></div>
</div>

<div class="tab" id="tab3">
  <h3 style="margin-top: 0px;">Animation looping the major and minor TSs</h3>

  <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden;">
    <video style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;" controls loop>
      <source src="/assets/electrophilic_s_amination/S_amination_annotated.mp4"></source>
      Your browser does not support the video tag.
    </video>
  </div>

</div>

<script src="/assets/tabs.js"></script>

---
{: data-content="footnotes"}

[^ref]: Read the full [paper](https://doi.org/10.1021/jacs.5c15618) on the J. Am. Chem. Soc. website
[^ACS_guidelines]: See [ACS guidelines for scholarly posting & sharing policies](https://pubs.acs.org/page/copyright/journals/posting_policies.html)