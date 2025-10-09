function toggleTab(panelId, tabIndex) {
  // Get the specific panel container
  let panel = document.getElementById(panelId);
  if (!panel) return;

  // Get buttons within this panel
  let buttons = panel.getElementsByClassName("tab-button");

  // Get all tabs in the document (like your original code)
  let allTabs = document.getElementsByClassName("tab");
  
  // Find the starting index for this panel's tabs
  // We'll use a data attribute or calculate based on panel order
  let panels = document.querySelectorAll(".tab-buttons");
  let startIndex = 0;
  
  for (let i = 0; i < panels.length; i++) {
    if (panels[i].id === panelId) {
      break;
    }
    startIndex += panels[i].getElementsByClassName("tab-button").length;
  }

  // Hide all tabs for this specific panel
  let tabCount = buttons.length;
  for (let i = 0; i < tabCount; i++) {
    if (allTabs[startIndex + i]) {
      allTabs[startIndex + i].style.display = "none";
    }
  }

  // Show the selected tab
  if (allTabs[startIndex + tabIndex]) {
    allTabs[startIndex + tabIndex].style.display = "block";
  }

  // Remove active class from all buttons in this panel
  for (let i = 0; i < buttons.length; i++) {
    buttons[i].classList.remove("active-tab-button");
  }

  // Add active class to the clicked button
  if (buttons[tabIndex]) {
    buttons[tabIndex].classList.add("active-tab-button");
  }
}

// Initialize all tab panels immediately
function initializeTabPanels() {
  let panels = document.getElementsByClassName("tab-buttons");
  for (let panel of panels) {
    toggleTab(panel.id, 0);
  }
}

// Run initialization immediately (like your original code)
initializeTabPanels();