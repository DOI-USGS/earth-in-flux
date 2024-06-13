<template>
    The route for this page is {{ pageRoute }}
    <ReferencesSection />
    <AuthorshipSection />
</template>

<script setup>
    import { useRoute } from 'vue-router'

    import ReferencesSection from '@/components/ReferencesSection.vue';
    import AuthorshipSection from '@/components/AuthorshipSection.vue';

    // global variables
    const route = useRoute()
    const pageRoute = route.params.pageRoute
    const pageContent = route.params.pageContent
</script>

<style scoped lang="scss">
</style>