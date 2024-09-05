<template>
  <section id="visualization-container">
    <div class="page-content">
      <div class="text-container" :class="{ mobile: mobileView}">
        <div v-if="!projectPage">
          <div class="typewriter">
            <h1 class="title">
              {{ text.pageTitle }}
            </h1>
          </div>
          <h2 class="subtitle">
            {{  text.pageSubTitle }}
          </h2>
        </div>
        <h2 v-if="projectPage" class="subtitle">
          {{ projectText.title }} project
        </h2>
      </div>
      <ChartGrid v-if="!projectPage" :view="currentView"/>
      <!---VizSection-->
      <VizSection
          v-if="projectPage"
          :figures="false"
          :fig-caption="false"
      >
          <!-- HEADING -->
          <template #heading>
            <h3>About the {{ projectText.title }} research</h3>
          </template>
          <template #aboveExplanation>
              <p v-html="projectText.motivation" />
          </template>
      </VizSection>
      <VizSection
          v-if="projectPage"
          :figures="true"
          :fig-caption="false"
      >
          <!-- HEADING -->
          <template #heading>
            <h3>{{ projectText.title }} visualizations</h3>
          </template>
          <template #figures>
            <ChartGrid :view="currentView"/>
          </template>
      </VizSection>
      <VizSection
          v-if="projectPage"
          :figures="true"
          :fig-caption="false"
      >
          <!-- HEADING -->
          <template #heading>
            <h3 v-if="projectText.teamData">Meet the team</h3>
          </template>
          <template #aboveExplanation>              
              <p v-if="projectText.teamData && projectText.teamText" v-html="projectText.teamText" />
          </template>
          <template #figures>
            <AboutTheTeam v-if="projectText.teamData" :key="projectRoute" :data="projectText.teamData"/>
          </template>
      </VizSection>
      <ReferencesSection v-if="projectReferences" title="Project resources" titleLevel="3" :references="projectReferences"/>
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
  import AboutTheTeam from '@/components/AboutTheTeam.vue';
  import ReferencesSection from '@/components/ReferencesSection.vue';
  import references from "@/assets/text/references";
  import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";

  // global variables
  const route = useRoute()
  const mobileView = isMobile;
  const projectRoute = ref(route.params.projectRoute)
  const projectKey = ref(projectRoute.value ? `${projectRoute.value.replace(/-/g, '')}` : null)
  const currentView = computed(() => {
    return projectRoute.value ? projectRoute.value : 'all'
  });
  const projectPage = computed(() => {
    return projectRoute.value ? true : false
  });
  const projectText = computed(() => {
    return projectRoute.value ? text.projects[projectKey.value] : null
  });
  const projectReferences = computed(() => {
    return projectRoute.value ? references[projectKey.value] : null
  })
  const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;

  //watches router params for changes
  watch(route, () => {
    projectRoute.value = route.params.projectRoute
    projectKey.value = projectRoute.value ? `${projectRoute.value.replace(/-/g, '')}` : null
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