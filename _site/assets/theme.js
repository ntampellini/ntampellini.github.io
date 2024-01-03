// Adapted from CodingBobby's theme

let theme = localStorage.getItem('theme')
let lightModeHTML = '☀'
let darkModeHTML = '☾'

function themeButton(theme) {
   let button = document.getElementById("theme-toggle")

   if (theme === 'light') {
      button.innerHTML = lightModeHTML
   } else if (theme === 'dark') {
      button.innerHTML = darkModeHTML
   }
}

function themeAttribute(theme) {
   if (theme === 'light') {
      document.documentElement.setAttribute('data-theme', 'light')
   } else if (theme === 'dark') {
      document.documentElement.setAttribute('data-theme', 'dark')
   }
}

function setTheme(theme) {
   if (theme === 'light') {
      themeAttribute('light')
      localStorage.setItem('theme', 'light')
      themeButton('dark')
   } else if (theme === 'dark') {
      themeAttribute('dark')
      localStorage.setItem('theme', 'dark')
      themeButton('light')
   }
}


function modeSwitcher() {

   let theme = localStorage.getItem('theme')

   if (theme === "dark") {
      setTheme('light')
   } else if (theme === "light") {
      setTheme('dark')
   }
}

// overly specific here to prevent usage of false storage contents
if (theme === "dark") {
   setTheme('dark')
} else {
   setTheme('light')
}
