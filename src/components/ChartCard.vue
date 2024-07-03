<template>
    <div class="chart">
        <img :src="src" :alt="alt"/>
    <div class="overlay">
        <p>{{ description }}</p>
    </div>
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
        description: {
            type: String,
            default: ''
        }
    })

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

    .overlay {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background-color: rgba(0, 0, 0, 0.5);
        color: white;
        opacity: 0;
        display: flex;
        align-items: center;
        justify-content: center;
    }

    .chart:hover .overlay {
        opacity: 1;
    }

    .overlay p {
        margin: 0;
        padding: 10px;
        text-align: center;
    }
</style>