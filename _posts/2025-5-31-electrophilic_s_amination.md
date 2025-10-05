---
layout: post
title: "[WIP] Catalytic asymmetric amination of sulfenamides"
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

This blog post contains interactive visualization for the "Catalytic asymmetric amination of sulfenamides to sulfinamidines enabled by first-principles catalyst design" <i>(manuscript in preparation)</i>.

<div>
  <img style='max-width: 95vw; width: 1000px; object-fit: contain' src="/assets/electrophilic_s_amination/reaction.jpg"/>
  <div align="center"><i>The electrophilic amination reaction of this work.</i></div>
</div>

### Major product TS (pro-<i>S</i>)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/electrophilic_s_amination/major.html"
    frameBorder="0"
  ></iframe>
</div>
<div align="center"><i>Substrate and aminating reagent highlighted in orange.</i></div>

### Minor product TS (pro-<i>R</i>)

<div class="iframe-parent">
  <iframe
    class="iframe"
    src="/assets/electrophilic_s_amination/minor.html"
    frameBorder="0"
  ></iframe>
</div>
<div align="center"><i>Substrate and aminating reagent highlighted in orange.</i></div>
