import markdownItKatex from 'markdown-it-katex';

export default {
    title: 'Meson Structure',
    description: 'Documentation for meson-structure analysis',
    base: '/meson-structure/',
    themeConfig: {
        sidebar: [
            { text: 'MC Variables', link: '/mc-variables' },
            { text: 'Data', link: '/data' },
            { text: 'Plots', link: '/plots' },
            { text: 'Analysis', link: '/analysis' }
        ]
    },
    head: [
        ['link', { rel: 'stylesheet', href: 'https://cdn.jsdelivr.net/npm/katex@0.16.0/dist/katex.min.css' }]
    ],
    markdown: {
        config: (md) => {
            md.use(markdownItKatex);
        }
    }
};