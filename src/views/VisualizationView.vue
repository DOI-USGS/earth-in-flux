<template>
  <section id="visualization-container">
    <div class="title-container" :class="{ mobile: mobileView}">
      <h1 class = 'title'>
        {{ text.pageTitle }}
      </h1>
    </div>
    <ChartGrid :view="currentView"/>
    <!---VizSection-->
    <VizSection
        v-if="projectPage"
        :figures="false"
        :fig-caption="false"
    >
        <!-- HEADING -->
        <template #heading>
            <h2>
              About the {{ projectBlurbText.title }} project
            </h2>
        </template>
        <template #aboveExplanation>
            <p v-html="projectBlurbText.blurb" />
        </template>
    </VizSection>
    <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
  </section>
</template>

<script setup>
  import { useRoute } from 'vue-router'
  import { isMobile } from 'mobile-device-detect';
  import text from "@/assets/text/text.js";
  import ChartGrid from '@/components/ChartGrid.vue';
  import VizSection from '@/components/VizSection.vue';
  import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

  // global variables
  const route = useRoute()
  const mobileView = isMobile;
  const projectRoute = route.params.projectRoute
  const currentView = projectRoute ? projectRoute : 'all'
  const projectPage = projectRoute ? true : false
  const projectBlurbText = projectRoute ? text.projects[`${projectRoute.replace(/-/g, '')}`] : null
  const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
</script>

<style scoped>
</style>