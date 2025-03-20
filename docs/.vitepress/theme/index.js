// docs/.vitepress/theme/index.js
import { h } from 'vue'
import DefaultTheme from 'vitepress/theme'
import './custom.css'

export default {
    extends: DefaultTheme,
    Layout: () => {
        return h(DefaultTheme.Layout, null, {
            // You can add custom layout slots here if needed
        })
    },
    enhanceApp({ app, router, siteData }) {
        // Register custom global components or other app enhancements
    }
}