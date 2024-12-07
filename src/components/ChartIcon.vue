<template>
    <div class="chart" :class="{ 'active': active }">
        <RouterLink :to="vizRoute ? { name: 'SubPage', params: { projectRoute, vizRoute } } : { name: 'ProjectPage', params: { projectRoute } }">
            <img v-if="vizRoute" :src="src" :alt="alt" />
            <div v-if="!vizRoute" id="overlay">
                <p class='overlay-text'>Project info</p>
            </div>
        </RouterLink>
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
        active: {
            type: Boolean,
            default: false
        },
        vizRoute: {
            type: String,
            default: ''
        },
        projectRoute: {
            type: String,
            default: ''
        },
    })
</script>

<style scoped lang="scss">
    $card-height: 90px;
	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: $card-height;
		height: $card-height;
        border-color: var(--faded-usgs-blue);
        border-width: 0px;
        border-style: solid;
        border-radius: 1.5px;
        box-shadow: 0px 0px 8px rgba(39,44,49,.07), 1px 4px 4px rgba(39,44,49,.04);
        opacity: 0.5;
	}
    .chart img {
        position: absolute;
        width: 100%;
        top: 0;
        background-color: white;
        border-radius: 1.5px;
    }
    
    .active {
        border-color: var(--faded-usgs-blue);
        border-width: 3px;
        border-style: solid;
        border-radius: 1.5px;
        opacity: 1;
        transform: scale(1.1);
    }

    .chart:hover {
        opacity: 1;
        transform: scale(1.1);
    }
    #overlay {
        height: $card-height;
        background-color: var(--faded-usgs-blue);
        border-radius: 3px;
    }
    .overlay-text {
        color: var(--usgs-blue);
        text-align: center;
        font-weight: 300;
        position: absolute;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        padding: 0;
    }
</style>