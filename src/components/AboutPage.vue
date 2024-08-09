<template>
        <div class="page-content">
            <!---VizSection-->
            <VizSection
                :figures="false"
                :fig-caption="false"
            >
                <!-- HEADING -->
                <template #heading>
                    <h2>
                        About the Climate Chart Gallery project
                    </h2>
                </template>
                <template #aboveExplanation>
                    <p v-if="text.collaborationDescription" v-html="text.collaborationDescription" />
                </template>
            </VizSection>
            <!---VizSection-->
            <!---VizSection-->
            <VizSection
                :figures="true"
                :fig-caption="false"
            >
                <!-- HEADING -->
                <template #heading>
                    <h3 v-if="vizlabText.teamData">The USGS Vizlab team</h3>
                </template>
                <template #aboveExplanation>
                    <p v-if="vizlabText.teamData && vizlabText.teamText" v-html="vizlabText.teamText" />
                </template>
                <template #figures>
                    <AboutTheTeam v-if="vizlabText.teamData" :data="vizlabText.teamData"/>
                </template>
            </VizSection>
            <!---VizSection-->
            <VizSection
                :figures="false"
                :fig-caption="false"
            >
                <!-- HEADING -->
                <template #heading>
                    <h3>
                        USGS climate projects
                    </h3>
                </template>
                <template #aboveExplanation>
                    <p>The current USGS climate projects highlighted in this page are the
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
        <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
</template>

<script setup>
    import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";
    import VizSection from '@/components/VizSection.vue';
    import AboutTheTeam from '@/components/AboutTheTeam.vue'
    import text from "@/assets/text/text.js";
    
    // global variables
    const vizlabText = text.vizlab;
    const projectKeys = Object.keys(text.projects).sort(); // Get alphabetical list of project keys
    const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
</script>

<style scoped lang="scss">
</style>