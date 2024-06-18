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
        <ReferencesSection :references="vizReferences"/>
        <AuthorshipSection />
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

    // global variables
    const route = useRoute()
    const vizRoute = route.params.vizRoute
    const projectRoute = route.params.projectRoute
    const mobileView = isMobile;
    const chartContent = ChartGrid.chartGridItems;
	const filteredChartContent = chartContent.filter(d => d.vizRoute === vizRoute)[0]
    const vizReferences = references[`${filteredChartContent.vizKey}`]

    const VizComponent = defineAsyncComponent(() =>
        import(`./${filteredChartContent.vizKey}Viz.vue`)
    )
</script>

<style scoped lang="scss">
</style>