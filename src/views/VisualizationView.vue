<template>
  <section id="visualization-container">
    <div class="page-content">
      <div class="text-container" :class="{ mobile: mobileView}">
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
    </div>
    <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
  </section>
</template>

<script setup>
  import { useRoute } from 'vue-router';
  import { computed, ref, watch } from 'vue';
  import { isMobile } from 'mobile-device-detect';
  import text from "@/assets/text/text.js";
  import ChartGrid from '@/components/ChartGrid.vue';
  import VizSection from '@/components/VizSection.vue';
  import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

  // global variables
  const route = useRoute()
  const mobileView = isMobile;
  const projectRoute = ref(route.params.projectRoute)
  const currentView = computed(() => {
      return projectRoute.value ? projectRoute.value : 'all'
  });
  const projectPage = computed(() => {
      return projectRoute.value ? true : false
  });
  const projectBlurbText = computed(() => {
      return projectRoute.value ? text.projects[`${projectRoute.value.replace(/-/g, '')}`] : null
  });
  
  const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;

  //watches router params for changes
  watch(route, () => {
    projectRoute.value = route.params.projectRoute
  })

</script>

<style scoped>
</style>