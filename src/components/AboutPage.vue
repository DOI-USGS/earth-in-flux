<template>
        <div class="page-content">
            <div
                class="text-container" 
                :class="{ mobile: mobileView}"
            >
                <p>Short blurb about <a href="https://labs.waterdata.usgs.gov/visualizations/index.html#/" target="_blank">VizLab</a> team / about the project</p>
            </div>
        </div>
        <PreFooterCodeLinks :gitHubRepositoryLink="gitHubRepositoryLink"/>
</template>

<script setup>
    import PreFooterCodeLinks from "@/components/PreFooterCodeLinks.vue";
    import { isMobile } from 'mobile-device-detect';
    
    const gitHubRepositoryLink = import.meta.env.VITE_APP_GITHUB_REPOSITORY_LINK;
    const mobileView = isMobile;
</script>

<style scoped lang="scss">
</style>