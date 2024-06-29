---
layout: post
title: Haiku Extractor
---
<style>

    .myButton {
        background-color:#ffd9ce;
        border-radius:17px;
        display:inline-block;
        cursor: pointer;
        color:#000000;
        padding: 7px 12px;
        text-decoration: none;
        border-width: 0px;
    }

    .myButton:hover {
        background-color:#e86f87;
        color:#000000;
    }

    .myButton:active {
        position: relative;
        top:1px;
    }

    .bottomGap {
        padding-bottom: 2vh;
    }

    textarea {
        color: var(--main-text-color);
        background-color: var(--quote-background);
        max-width: min(70ch, 80vw);
        min-width: min(70ch, 80vw);
        width: min(70ch, 80vw);
        height: 15vw;
        padding: 0.25rem;
        font-family: Fira Code, monospace;
        border-width: 3px;
        text-align: center;
    }

    .haiku {
        text-align: center;
        font-size: 2rem;
    }

    .haiku2 {
        text-align: center;
        /* text-color: var(--quote-color); */
        opacity: 0.33;
        font-style: italic;
        font-size: 1.5rem;
    }

</style>
        
Paste a block of text in the box below to extract random Haiku pieces from it.


<!-- <html>
    <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width,initial-scale=1" />
        <link rel="stylesheet" href="https://pyscript.net/releases/2024.1.1/core.css" />
        <script type="module" src="https://pyscript.net/releases/2024.1.1/core.js"></script>
    </head>
    
    <body> -->

<script type="module" src="https://pyscript.net/releases/2024.1.1/core.js"></script>

<div align="center" class="bottomGap">
    <textarea wrap="soft" name="text" id="text" placeholder="Paste a block of text..."></textarea>
</div>

<div align="center" class="bottomGap">
    <label class="container" align="center">Contiguous verses
    <input id="checkbox" type="checkbox" checked="checked">
    <span class="checkmark"></span>
    </label>
</div>

<div align="center">
    <button py-click="main" class="myButton">Get haiku</button>
    <button py-click="scrape_chemrxiv" class="myButton">Scrape ChemRXiv</button>
</div>

<br>
<div id="paper_credits" class="haiku2"></div>

<br>
<div id="output" class="haiku"></div>
<br>

<div id="output_extra" class="haiku2"></div>

<script type="py" src="../assets/haiku/haiku.py" config="../assets/haiku/pyscript.json"></script>