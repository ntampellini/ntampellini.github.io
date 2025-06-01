---
layout: post
title: "Enantiocontrolled Cyclization of Inherently Chiral 7- and 8-Membered Rings"
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
</style>

<!-- CONTENT -->

> Click and drag to rotate, Ctrl+click/Mouse3 to move, scroll/right-click to zoom.

This blog post contains interactive visualization for the "Enantiocontrolled Cyclization to Form Chiral 7- and 8‑Membered Rings Unified by the Same Catalyst Operating with Different Mechanisms".[^JACS]

<div>
  <img style='max-width: 95vw; width: 700px; object-fit: contain' src="/assets/enantiocontrolled_cyclization/reaction.jpg"/>
  <div align="center"><i>Graphical TOC of the paper.</i></div>
</div>

In the following sections you can find the best conformers for the transition states leading to the major enantiomers of the representative 7- and 8-membered ring substrates. Surprisingly, even though they are formed via the *same* catalyst, they feature *opposite* absolute configurations! You can rotate the structures and toggle the visualization of the dipole analysis, which shows how the local components of the catalyst and substrate appear to be oriented in such a way to reduce the magnitude of the overall dipole. This effect, along with hydrogen bonding and dipersion interactions, dictates the relative arrangement of the transition state components and ultimately enables selectivity. The dipole decomposition was performed using atomic partial charges computed with the MBIS method[^MBIS] through ORCA[^ORCA] on DFT-optimized structures. Interactive visuals were created with Python and Javascript (3DMol.js[^3DMol]).

### S<sub>N</sub>Ar, major enantiomer TS (<i>M</i> configuration)

<div class="iframe-parent">
  <iframe
    id="7mem"
    class="iframe"
    src="/assets/enantiocontrolled_cyclization/7_major_dipole.html"
    frameBorder="0"
  ></iframe>
</div>

<div align="center" class="bottomGap">
    <label class="container" align="center"> Toggle Dipole analysis
    <input id="7mem_button" type="checkbox" checked="checked">
    <span class="checkmark"></span>
    </label>
</div>

### S<sub>N</sub>2, major enantiomer TS (<i>P</i> configuration)

<div class="iframe-parent">
  <iframe
    id="8mem"
    class="iframe"
    src="/assets/enantiocontrolled_cyclization/8_major_dipole.html"
    frameBorder="0"
  ></iframe>
</div>

<div align="center" class="bottomGap">
    <label class="container" align="center"> Toggle Dipole analysis
    <input id="8mem_button" type="checkbox" checked="checked">
    <span class="checkmark"></span>
    </label>
</div>

<script>

    var iframes = ["7mem", "8mem"];
    iframes.forEach(function(iframeName) {

      const checkbox = document.getElementById(iframeName+'_button');
      const iframe = document.getElementById(iframeName);
      
      checkbox.addEventListener('change', function() {
          // Send message to iframe
          iframe.contentWindow.postMessage({
              type: 'toggleButton',
              show: this.checked
          });
      });
      
      // Optional: Wait for iframe to load before sending initial state
      iframe.addEventListener('load', function() {
          setTimeout(() => {
              iframe.contentWindow.postMessage({
                  type: 'toggleButton',
                  show: checkbox.checked
              });
          }, 1000); // Give 3DMol time to initialize
      });
    });

</script>

---
{: data-content="footnotes"}

[^JACS]: The final version of this work is now published in [JACS](https://pubs.acs.org/doi/10.1021/jacs.4c17080).
[^MBIS]: See the original reference for MBIS charges [here](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00456).
[^ORCA]: ORCA is a wonderful quantum chemistry program, see the [ORCA forum](https://orcaforum.kofo.mpg.de/app.php/portal).
[^3DMol]: 3DMol.js is a powerful library for molecular visualization, see the [homepage](https://3dmol.csb.pitt.edu/) of the project.