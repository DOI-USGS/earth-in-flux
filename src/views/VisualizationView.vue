<template>
  <section id="visualization-container">
    <div class="page-content">
      <div class="text-container" :class="{ mobile: mobileView}">
        <div v-if="!projectPage">
          <div class="typewriter">
            <h1 class="title">
              {{ pageText.pageTitle }}
            </h1>
          </div>
          <h2 class="subtitle">
            {{  pageText.pageSubTitle }}
          </h2>
        </div>
        <h2 v-if="projectPage" class="subtitle">
          The {{ pageText.title }} project
        </h2>
      </div>
      <ChartGrid v-if="!projectPage" :view="currentView"/>
      <div v-if="!projectPage">
        <!---About the project-->
        <VizSection
          :figures="false"
          :fig-caption="false"
        >
          <template #heading>
            <hr>
            <h2>
              {{ pageText.collaborationHeading }}
            </h2>
          </template>
          <template #aboveExplanation>
            <p v-if="pageText.collaborationDescription" v-html="pageText.collaborationDescription" />
          </template>
        </VizSection>
        <!---About Vizlab-->
        <VizSection
          :figures="true"
          :fig-caption="false"
        >
          <template #heading>
            <h3 v-if="pageText.teamData"> {{ pageText.teamHeading }} </h3>
          </template>
          <template #aboveExplanation>
            <p v-if="pageText.teamData && pageText.teamText" v-html="pageText.teamText" />
          </template>
          <template #figures>
            <AboutTheTeam v-if="pageText.teamData" :data="pageText.teamData"/>
          </template>
        </VizSection>
        <!---Highlighted Projects-->
        <VizSection
          :figures="false"
          :fig-caption="false"
        >
          <template #heading>
            <h3>
              {{ pageText.projectsHeading }}
            </h3>
          </template>
          <template #aboveExplanation>
            <p> {{ pageText.projectsLeadIn }}
            <span v-for="projectKey in projectKeys" :key="projectKey">
              <RouterLink :key="projectKey" :to="{ path: `/${text.projects[projectKey].title.replace(/\s+/g, '-').toLowerCase()}` }"> {{ text.projects[projectKey].title }}</RouterLink>
                <span v-if="projectKeys.indexOf(projectKey) != projectKeys.length - 1 && projectKeys.length > 2">, </span>
                <span v-if="projectKeys.indexOf(projectKey) == projectKeys.length - 2"> and </span>
            </span>
            projects.
            </p>
          </template>
        </VizSection>
      </div>
      <!---Project page-->
      <div v-if="projectPage">
        <VizSection
          :figures="false"
          :fig-caption="false"
        >
          <template #heading>
            <h3>About the {{ pageText.title }} research</h3>
          </template>
          <template #aboveExplanation>
            <div v-for="paragraph, item in pageText.motivation" :key="item">
              <p v-html="paragraph" />
            </div>
          </template>
        </VizSection>
        <VizSection
          :figures="true"
          :fig-caption="false"
        >
          <template #heading>
            <h3>{{ pageText.title }} visualizations</h3>
          </template>
          <template #figures>
            <ChartGrid :view="currentView"/>
          </template>
        </VizSection>
        <VizSection
          :figures="true"
          :fig-caption="false"
        >
          <template #heading>
            <h3 v-if="pageText.teamData">Meet the team</h3>
          </template>
          <template #aboveExplanation>              
            <p v-if="pageText.teamData && pageText.teamText" v-html="pageText.teamText" />
          </template>
          <template #figures>
            <AboutTheTeam v-if="pageText.teamData" :key="projectRoute" :data="pageText.teamData"/>
          </template>
        </VizSection>
      </div>
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
  const projectKeys = Object.keys(text.projects).sort(); // Get alphabetical list of project keys
  const currentView = computed(() => {
    return projectRoute.value ? projectRoute.value : 'all'
  });
  const projectPage = computed(() => {
    return projectRoute.value ? true : false
  });
  const pageText = computed(() => {
    return projectRoute.value ? text.projects[projectKey.value] : text.landingPage
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