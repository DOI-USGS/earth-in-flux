<template>
    <section 
        :id="vizRoute"
    >
        <div
            class="text-container" 
            :class="{ mobile: mobileView}"
        >
            <h1 class = 'title'>
                {{ filteredChartContent.title }}
            </h1>
            The route for this page is {{ vizRoute }}. This viz is associated with the {{ projectRoute }} project.
        </div>
        <VizComponent/>
        <ReferencesSection v-if="vizReferences" :references="vizReferences"/>
        <AuthorshipSection v-if="vizAuthors" :authors="vizAuthors"/>
    </section>
</template>

<script setup>
    import { useRoute } from 'vue-router'
    import { isMobile } from 'mobile-device-detect';
    import { defineAsyncComponent } from 'vue'

    import ChartGrid from '@/assets/content/ChartGrid.js';
    import ReferencesSection from '@/components/ReferencesSection.vue';
    import references from "@/assets/text/references";
    import AuthorshipSection from '@/components/AuthorshipSection.vue';
    import authors from "@/assets/text/authors";

    // global variables
    const route = useRoute()
    const vizRoute = route.params.vizRoute
    const projectRoute = route.params.projectRoute
    const mobileView = isMobile;
    const chartContent = ChartGrid.chartGridItems;
	const filteredChartContent = chartContent.filter(d => d.vizRoute === vizRoute)[0]
    const vizReferences = references[`${filteredChartContent.vizKey}`]
    const vizAuthors = authors[`${filteredChartContent.vizKey}`]

    const VizComponent = defineAsyncComponent(() =>
        //Use convention of adding Viz to component names to ensure are two words
        //Plus, also allows us to specify a filename pattern here, per:
        //https://github.com/rollup/plugins/tree/master/packages/dynamic-import-vars#imports-to-your-own-directory-must-specify-a-filename-pattern
        import(`./${filteredChartContent.vizKey}Viz.vue`) //
    )
</script>

<style scoped lang="scss">
</style>