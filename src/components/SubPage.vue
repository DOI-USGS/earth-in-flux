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
            </div>
            <VizComponent :id="`${vizRoute}-viz`" :text="vizText" :key="$route.fullPath"/>
            <ChartIconGrid :projectRoute="projectRoute" :vizRoute="vizRoute" :key="vizRoute"/>
            <ReferencesSection v-if="vizReferences" title="References" titleLevel="3" :references="vizReferences"/>
            <AuthorshipSection v-if="vizAuthors" title="" titleLevel="3" :authors="vizAuthors"/>
        </div>
        <PreFooterCodeLinks :gitHubRepositoryLink="vizGitHubRepositoryLink"/>
    </section>
</template>

<script setup>
    import { useRoute } from 'vue-router'
    import { isMobile } from 'mobile-device-detect';
    import { computed, defineAsyncComponent, ref, shallowRef, watch } from 'vue'
    import ChartIconGrid from '@/components/ChartIconGrid.vue';
    import ChartGridContent from '@/assets/content/ChartGrid.js';
    import text from "@/assets/text/text.js";
    import ReferencesSection from '@/components/ReferencesSection.vue';
    import references from "@/assets/text/references";
    import AuthorshipSection from '@/components/AuthorshipSection.vue';
    import authors from "@/assets/text/authors";    
    import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

    // global variables
    const route = useRoute()
    const vizRoute = ref(route.params.vizRoute);
    const projectRoute = ref(route.params.projectRoute);
    const mobileView = isMobile;
    const VizComponent = shallowRef(
        defineAsyncComponent(() =>
            //Use convention of adding Viz at end of component names to ensure are two words
            //Plus, also allows us to specify a filename pattern here, per:
            //https://github.com/rollup/plugins/tree/master/packages/dynamic-import-vars#imports-to-your-own-directory-must-specify-a-filename-pattern
            import(`./${vizKey.value}Viz.vue`) //
        )
    );
    const chartContent = ChartGridContent.chartGridItems;
    const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
    
    // set up computed properties based on ref values
    const filteredChartContent = computed(() => {
        return chartContent.filter(d => d.vizRoute === vizRoute.value)[0];
    }); 

    const vizKey = computed(() => {
        return filteredChartContent.value.vizKey;
    }); 

    const vizText = computed(() => {
        return text.visualizations[`${filteredChartContent.value.vizKey}`]
    }); 

    const vizReferences = computed(() => {
        return references[`${filteredChartContent.value.vizKey}`]
    }); 

    const vizAuthors = computed(() => {
        return authors[`${filteredChartContent.value.vizKey}`]
    }); 

    const vizGitHubRepositoryLink = computed(() => {
        return `${gitHubRepositoryLink}/blob/main/src/components/${vizKey.value}Viz.vue` //Use convention of adding Viz at end of component names to ensure are two words
    }); 

    //watches router params for changes
    watch(route, () => {
        // update project route
        projectRoute.value = route.params.projectRoute;
        // update viz route
        vizRoute.value = route.params.vizRoute;
        // update dynamic component
        VizComponent.value = defineAsyncComponent(() =>
            //Use convention of adding Viz at end of component names to ensure are two words
            //Plus, also allows us to specify a filename pattern here, per:
            //https://github.com/rollup/plugins/tree/master/packages/dynamic-import-vars#imports-to-your-own-directory-must-specify-a-filename-pattern
            import(`./${vizKey.value}Viz.vue`) //
        )
    })
</script>

<style scoped lang="scss">
</style>