<template>
  <section id="visualization-container">
    <div class="page-content">
      <div class="text-container" :class="{ mobile: mobileView}">
        <div class="typewriter">
          <h1 class="title">
            {{ text.pageTitle }}
          </h1>
        </div>
        <h2 v-if="projectPage" class="subtitle">
          {{ projectBlurbText.title }} project visualizations
        </h2>
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
  /* https://css-tricks.com/snippets/css/typewriter-effect/ */
  .typewriter h1 {
    overflow: hidden; /* Ensures the content is not revealed until the animation */
    border-right: .15em solid var(--usgs-blue); /* The typwriter cursor */
    white-space: nowrap; /* Keeps the content on a single line */
    margin: 0 auto; /* Gives that scrolling effect as the typing happens */
    letter-spacing: .1em; /* Adjust as needed */
    animation: 
      typing 3.5s steps(40, end),
      blink-caret .9s step-end infinite;
  }

  /* The typing effect */
  @keyframes typing {
    from { width: 0 }
    to { width: 100% }
  }

  /* The typewriter cursor effect */
  @keyframes blink-caret {
    from, to { border-color: transparent }
    50% { border-color: var(--usgs-blue); }
  }
</style>