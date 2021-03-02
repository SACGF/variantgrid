/**
 *
 * Changed dlawrence 21/Nov/2014 - Use activeClass to move groups around not passed in groups
 *
 * Based on: JQuery MultiDraggable Plugin
 *
 * Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php)
 *
 * Written by Sudheer Someshwara <sudheer.someshwara@gmail.com>
 *
 * MultiDraggable is a jQuery plugin which extends jQuery UI Draggable to add multi drag and live functionality.
 *
 * Used version edited by j.verdurmen for jsplumb test - http://jsfiddle.net/sporritt/ZQcEc/25/
**/
(function ($, undefined) {
    $.fn.multiDraggable = function (opts) {
        var initLeftOffset = []
        , initTopOffset = [];
        return this.each(function () {
			function getItems(self) {
				if ($(self).hasClass(opts.activeClass)) {
					return $("." + opts.activeClass);
				}
				return $();
			}

            $(this).data("init", true).draggable(opts, {
                start: function (event, ui) {
                    var item = $(this);
                    var pos = item.position();
				    var items = getItems(this);
                    //console.log("start",ui,event)
					if (items.length) {
                        $.each(items, function (key, value) {
                            var elemPos = $(value).position();
                            initLeftOffset[key] = elemPos.left - pos.left;
                            initTopOffset[key] = elemPos.top - pos.top;
                            opts.startAll ? opts.startAll.call(this, event, ui) : {};
                        });
                        jsPlumb.repaint(items);
                    }
                    opts.startNative ? opts.startNative.call(this, event, ui) : {};
                    //  items.trigger("start");
                    jsPlumb.repaint(item);
                },
				drag: function (event, ui) {
				    var item = $(this);
				    var pos = ui.offset;
				    var items = getItems(this);
					if (items.length) {
					    $.each(items, function (key, value) {
							var oPos = {
							    left: pos.left + initLeftOffset[key],
							    top: pos.top + initTopOffset[key]
							}, oEl = $(value);
						    
							oEl.offset(oPos);
							jsPlumb.repaint(oEl, oPos);
					    });
					}
				    // repaint the dragging item, passing in position                            
				    opts.dragNative ? opts.dragNative.call(this, event, ui) : {};                            
				    jsPlumb.repaint(item, pos);
				},
                stop: function (event, ui) {
                    var item = $(this);
                    var pos = $(this).offset();
					var items = getItems(this);
             //         console.log("stop",ui,event)
					if (items.length) {
                        $.each(items, function (key, value) {
                            $(value).offset({
                                left: pos.left + initLeftOffset[key],
                                top: pos.top + initTopOffset[key]
                            });
                            opts.stopAll ? opts.stopAll.call(this, event, ui) : {};
                        });
                        jsPlumb.repaint(items);
                    }
                    opts.stopNative ? opts.stopNative.call(this, event, ui) : {};
                    //     items.trigger("stop");
                    //jsPlumb.repaint(item);
                },
            });
        });
    };
}(jQuery));
