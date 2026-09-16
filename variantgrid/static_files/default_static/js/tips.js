// Feature tips - server side is variantgrid/tips.py, markup is tips/tip_box.html.
// A box shuffles its tip list once, then the next-tip button walks through it, so a user sees the
// whole set before anything repeats.

function shuffleTips(tips) {
    const shuffled = tips.slice();
    for (let i = shuffled.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    return shuffled;
}

function turnOffTips() {
    $.ajax({
        type: 'POST',
        url: Urls.set_show_tips(),
        headers: {'X-CSRFToken': Cookies.get('csrftoken')}
    });
    const message = $("<div>", {class: "tip-box tips-off", text: "Tips off - switch them back on in User Settings."});
    const $boxes = $(".tip-box");
    $boxes.first().replaceWith(message);
    $(".tip-box").not(message).remove();
    message.delay(4000).fadeOut();
}

function setupTipBox(box) {
    const $box = $(box);
    if ($box.data("tipsSetup")) {
        return;  // boxes get re-scanned when the analysis grid overlay adds one
    }
    $box.data("tipsSetup", true);

    // Start from the tip the server already picked, so the button moves us on rather than back
    const shuffled = shuffleTips($box.data("tips") || []);
    let index = Math.max(shuffled.indexOf($box.find(".tip-text").text()), 0);

    $box.find(".tip-next").click(function() {
        if (shuffled.length < 2) {
            return;
        }
        index = (index + 1) % shuffled.length;
        $box.find(".tip-text").text(shuffled[index]);
    });

    $box.find(".tip-dismiss").click(turnOffTips);
}

function setupTipBoxes() {
    $(".tip-box").each(function() {
        setupTipBox(this);
    });
}

// Build a tip box for somewhere there's no template render - eg the analysis grid loading overlay,
// which gets its tips up front as ANALYSIS_TIPS. Returns null when there's nothing to show.
function buildTipBox(tips) {
    if (!tips || !tips.length) {
        return null;
    }
    const $box = $("<div>", {class: "tip-box"});
    $("<i>", {class: "fa-solid fa-lightbulb tip-icon"}).appendTo($box);
    $("<span>", {class: "tip-label", text: "Tip:"}).appendTo($box);
    $("<span>", {class: "tip-text", text: tips[Math.floor(Math.random() * tips.length)]}).appendTo($box);
    const $next = $("<a>", {class: "tip-next", href: "javascript:void(0)", title: "Show another tip"});
    $("<i>", {class: "fa-solid fa-arrows-rotate"}).appendTo($next);
    $next.appendTo($box);
    const $dismiss = $("<a>", {class: "tip-dismiss", href: "javascript:void(0)",
                               title: "Turn tips off - you can switch them back on in User Settings"});
    $("<i>", {class: "fa-solid fa-xmark"}).appendTo($dismiss);
    $dismiss.appendTo($box);

    $box.data("tips", tips);
    setupTipBox($box);
    return $box;
}

$(document).ready(function() {
    setupTipBoxes();
});
