<template>
  <section id="visualization-container">
    <div class="text-container" :class="{ mobile: mobileView}">
      <h1 class = 'title'>
        {{ text.pageTitle }}
      </h1>
    </div>
    <ChartGrid :view="currentView"/>
    <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
  </section>
</template>

<script setup>
  import { useRoute } from 'vue-router'
  import { isMobile } from 'mobile-device-detect';
  import text from "@/assets/text/text.js";
  import ChartGrid from '@/components/ChartGrid.vue';
  import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

  // global variables
  const route = useRoute()
  const mobileView = isMobile;
  const projectRoute = route.params.projectRoute
  const currentView = projectRoute ? projectRoute : 'all'
  const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
</script>

<style scoped>
</style>