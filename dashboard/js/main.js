/* main.js, register per-tab inits and boot. */
document.addEventListener("DOMContentLoaded", () => {
  Theme.init();
  TabManager.registerInit("overview", Overview.init);
  TabManager.registerInit("atlas", Atlas.init);
  TabManager.registerInit("scores", Scores.init);
  TabManager.registerInit("compare", Compare.init);
  TabManager.registerInit("within", Within.init);
  TabManager.registerInit("finder", Finder.init);
  TabManager.registerInit("about", About.init);
  TabManager.init();
});
