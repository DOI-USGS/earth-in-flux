<template>
    <div class="chart" v-if="show">
        <img :src="src" :alt="alt"/>
    </div>
</template>

<script setup>
    const props = defineProps({
        src: {
            type: String,
            default: ``
        },
        alt: {
            type: String,
            default: ``
        },
        project: {
            type: String,
            default: ``
        },
        view: {
            type: String,
            default: ``
        },
    })

    // global variables
    const projectRoute = props.project.replace(/\s+/g, '-').toLowerCase();
    let show = props.view == projectRoute ? true : false;
    show = props.view == 'all' ? true : show;
</script>

<style scoped lang="scss">
	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: 280px;
		height: 100%;
	}
    .chart:hover {
        transform: scale(1.1);
    }
</style>