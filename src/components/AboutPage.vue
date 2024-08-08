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
                    <AboutTheTeam v-if="vizlabText.teamData" :key="projectRoute" :data="vizlabText.teamData"/>
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
                    <p v-for="project in text.projects" :key="project">
                        <RouterLink :key="project" :to="{ path: `/${project.title.replace(/\s+/g, '-').toLowerCase()}` }"> {{ project.title }} </RouterLink>
                    </p>
                </template>
            </VizSection>
        </div>
        <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
</template>

<script setup>
    import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";
    // import { isMobile } from 'mobile-device-detect';
    import VizSection from '@/components/VizSection.vue';
    import AboutTheTeam from '@/components/AboutTheTeam.vue'
    import text from "@/assets/text/text.js";
    
    // global variables
    const vizlabText = text.vizlab
    const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
    // const mobileView = isMobile;
</script>

<style scoped lang="scss">
</style>