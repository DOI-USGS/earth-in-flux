<template>
    <header id="NavBar" class="header-nav"  role="banner">
        <div class="topnav">
            <RouterLink :to="{ path: '/' }">Home</RouterLink>
            <div class="dropdown">
                <a href="javascript:void(0)" class="dropdown">Projects</a>
                <div class="dropdown-content">
                    <RouterLink v-for="projectKey in projectKeys" :key="projectKey" :to="{ path: `/${text.projects[projectKey].title.replace(/\s+/g, '-').toLowerCase()}` }"> {{ text.projects[projectKey].title }} </RouterLink>
                </div>
            </div>
        </div>
    </header>
</template>

<script setup>
    import text from "@/assets/text/text.js";

    // Get alphabetical list of project keys
    const projectKeys = Object.keys(text.projects).sort()
</script>

<style scoped>

.header-nav {
    position: -webkit-sticky; /* For Safari */
    position: sticky;
    top: 0;
    z-index: 1000;
}

.topnav {
    background-color: #65778B;
    overflow: hidden;
    display: flex;
    justify-content: center;
}

.topnav a {
    color: #ffffff;
    text-align: center;
    padding: 14px 16px;
    text-decoration: none;
    font-weight: 400;
}

.topnav a:hover {
    background-color: var(--faded-usgs-blue);
    color: var(--usgs-blue);
}

.dropdown {
    float: left;
    overflow: hidden;
}

.dropdown .dropbtn {
    font-size: 16px;  
    border: none;
    outline: none;
    color: white;
    padding: 14px 16px;
    background-color: inherit;
    font-family: inherit;
    margin: 0;
}

.dropdown-content a:hover {
    background-color: var(--faded-usgs-blue);
    color: var(--usgs-blue);
}

.dropdown:hover .dropdown-content {
    display: block;
}

.dropdown-content {
    display: none;
    position: absolute;
    top: 48px;
    background-color: var(--color-background);
    min-width: 160px;
    box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
    z-index: 1;
}

.dropdown-content a {
    float: none;
    color: var(--usgs-blue);
    padding: 12px 16px;
    text-decoration: none;
    display: block;
    text-align: left;
}
</style>
