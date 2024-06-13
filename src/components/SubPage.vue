<template>
    The route for this page is {{ subPageRoute }}
    <ReferencesSection />
    <AuthorshipSection />
</template>

<script setup>
    import { useRoute } from 'vue-router'

    import ReferencesSection from '@/components/ReferencesSection.vue';
    import AuthorshipSection from '@/components/AuthorshipSection.vue';

    // global variables
    const route = useRoute()
    const subPageRoute = route.params.route
</script>

<style scoped lang="scss">
</style>