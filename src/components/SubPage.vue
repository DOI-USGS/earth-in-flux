<template>
    <section 
        :id="`${vizRoute}-subpage`"
    >
        <div class="page-content">
            <div
                class="text-container" 
                :class="{ mobile: mobileView}"
            >
                <h1 class = 'title'>
                    {{ filteredChartContent.title }}
                </h1>
                <p>This visualization is associated with the <RouterLink :to="`/${projectRoute}`"> {{ filteredChartContent.project }}</RouterLink> project.</p>
            </div>
            <VizComponent :id="`${vizRoute}-viz`" :text="vizText"/>
            <ReferencesSection v-if="vizReferences" title="References" titleLevel="2" :references="vizReferences"/>
            <AuthorshipSection v-if="vizAuthors" :authors="vizAuthors"/>
        </div>
        <PreFooterCodeLinks :gitHubRepositoryLink="vizGitHubRepositoryLink"/>
    </section>
</template>

<script setup>
    import { useRoute } from 'vue-router'
    import { isMobile } from 'mobile-device-detect';
    import { defineAsyncComponent } from 'vue'

    import ChartGrid from '@/assets/content/ChartGrid.js';
    import text from "@/assets/text/text.js";
    import ReferencesSection from '@/components/ReferencesSection.vue';
    import references from "@/assets/text/references";
    import AuthorshipSection from '@/components/AuthorshipSection.vue';
    import authors from "@/assets/text/authors";    
    import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

    // global variables
    const route = useRoute()
    const vizRoute = route.params.vizRoute
    const projectRoute = route.params.projectRoute
    const mobileView = isMobile;
    const chartContent = ChartGrid.chartGridItems;
	const filteredChartContent = chartContent.filter(d => d.vizRoute === vizRoute)[0]
    const vizKey = filteredChartContent.vizKey
    const vizText = text.visualizations[`${filteredChartContent.vizKey}`]
    const vizReferences = references[`${filteredChartContent.vizKey}`]
    const vizAuthors = authors[`${filteredChartContent.vizKey}`]
    const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
    const vizGitHubRepositoryLink = `${gitHubRepositoryLink}/blob/main/src/components/${vizKey}Viz.vue` //Use convention of adding Viz at end of component names to ensure are two words
    
    const VizComponent = defineAsyncComponent(() =>
        //Use convention of adding Viz at end of component names to ensure are two words
        //Plus, also allows us to specify a filename pattern here, per:
        //https://github.com/rollup/plugins/tree/master/packages/dynamic-import-vars#imports-to-your-own-directory-must-specify-a-filename-pattern
        import(`./${vizKey}Viz.vue`) //
    )
</script>

<style scoped lang="scss">
</style>