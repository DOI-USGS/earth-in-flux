<template>
    <section 
        :id="pageRoute"
    >
        <div
            class="text-container" 
            :class="{ mobile: mobileView}"
        >
            The route for this page is {{ pageRoute }}. This viz is associated with the {{ projectRoute }} project.
        </div>
        <ReferencesSection />
        <AuthorshipSection />
    </section>
</template>

<script setup>
    import { useRoute } from 'vue-router'
    import { isMobile } from 'mobile-device-detect';

    import ReferencesSection from '@/components/ReferencesSection.vue';
    import AuthorshipSection from '@/components/AuthorshipSection.vue';

    // global variables
    const route = useRoute()
    const pageRoute = route.params.pageRoute
    const projectRoute = route.params.projectRoute
    const mobileView = isMobile;
</script>

<style scoped lang="scss">
</style>