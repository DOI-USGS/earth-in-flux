<template>
    <div class="chart">
        <div class="takeaway">
            <p>{{ description }}</p>
        </div>
        <img :src="src" :alt="alt"/>
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
    $collapsed-height: 280px;
    $expanded-height: 400px;
    $height-difference: $expanded-height - $collapsed-height;
	.chart{
		display: flex;
		flex-direction: column-reverse;
		justify-content: flex-start;
		position: relative;
		width: $collapsed-height;
		height: $collapsed-height;
        -moz-transition: height 1s ease, border-width 1s ease, border-radius 1s ease;
        -webkit-transition: height 1s ease, border-width 1s ease, border-radius 1s ease;
        -o-transition: height 1s ease, border-width 1s ease, border-radius 1s ease;
        transition: height 1s ease, border-width 1s ease, border-radius 1s ease;
        background-color: var(--faded-usgs-blue);
        border-color: var(--faded-usgs-blue);
        border-width: 0px;
        border-style: solid;
        border-radius: 1.5px;
        box-shadow: 0px 0px 8px rgba(39,44,49,.07), 1px 4px 4px rgba(39,44,49,.04);
	}
    .chart:hover {
        height: $expanded-height;
        transform: scaleY(1);
        border-width: 4px;
        border-radius: 3px;
    }
    .expanded {
        height: $expanded-height;
        transform: scaleY(1);
        border-width: 4px;
        border-radius: 3px;
    }
    .takeaway {
        position: absolute;
        left: 0;
        right: 0;
        bottom: calc($height-difference / 2);
        line-height: 2rem;
        -webkit-transition: top 1s ease; 
        -moz-transition: top 1s ease; 
        -ms-transition: top 1s ease; 
        -o-transition: top 1s ease; 
        transition: top 1s ease;
    }
    .takeaway p {
        font-family: sans-serif; /* This is fallback font for old browsers */
        font-family: var(--title-font);
        color: var(--color-title-text);
        transform: translateY(50%);
        padding: 0.2rem 2rem 0.2rem 2rem;
        opacity: 0;
        transition: all 1s ease; 
    }
    .chart:hover .takeaway p {
        opacity: 1;
    }
    .expanded .takeaway p {
        opacity: 1;
    }
    .chart img {
        position: absolute;
        width: 100%;
        top: 0;
        background-color: white;
        border-radius: 1.5px;
        transition: all 1s ease;
    }
    .chart:hover img {
        border-radius: 3px;
    }
    .expanded img {
        border-radius: 3px;
    }
</style>