// docs/.vitepress/theme/index.js
import DefaultTheme from 'vitepress/theme'
import './custom.css'

export default {
    ...DefaultTheme,
    enhanceApp({ app }) {
        // register global components or other app enhancements
    }
}