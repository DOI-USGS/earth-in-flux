<template>
    <div class="chart">
        <img :src="getThumb(src)" :alt="alt"/>
    </div>
</template>

<script setup>
    defineProps({
        src: {
            type: String,
            default: ``
        },
        alt: {
            type: String,
            default: ``
        },
    })

    function getThumb(pic) {
        return 'https://labs.waterdata.usgs.gov/visualizations/thumbnails/'+pic
    }
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