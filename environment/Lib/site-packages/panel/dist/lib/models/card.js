import { Column, ColumnView } from "@bokehjs/models/layouts/column";
import * as DOM from "@bokehjs/core/dom";
import { classes, empty } from "@bokehjs/core/dom";
import { Column as ColumnLayout } from "@bokehjs/core/layout/grid";
import * as p from "@bokehjs/core/properties";
export class CollapseableColumnLayout extends ColumnLayout {
    constructor(items, collapsed = false) {
        super(items);
        this.collapsed = collapsed;
    }
    _measure_totals(row_heights, col_widths) {
        // Do not count height of collapsed height
        const heights = this.collapsed ? row_heights.slice(0, 1) : row_heights;
        const last = heights.length - 1;
        if (!this.collapsed)
            heights[last] = heights[last] + 5;
        return super._measure_totals(heights, col_widths);
    }
}
CollapseableColumnLayout.__name__ = "CollapseableColumnLayout";
export class CardView extends ColumnView {
    connect_signals() {
        super.connect_signals();
        this.connect(this.model.properties.collapsed.change, () => this._collapse());
        const { active_header_background, header_background, header_color } = this.model.properties;
        this.on_change([active_header_background, header_background, header_color], () => this.render());
    }
    _update_layout() {
        const items = this.child_views.map((child) => child.layout);
        this.layout = new CollapseableColumnLayout(items, this.model.collapsed);
        this.layout.rows = this.model.rows;
        this.layout.spacing = [this.model.spacing, 0];
        const sizing = this.box_sizing();
        if (this.model.collapsed)
            sizing.height = undefined;
        this.layout.set_sizing(sizing);
    }
    render() {
        empty(this.el);
        const { background, button_css_classes, header_color, header_tag, header_css_classes } = this.model;
        this.el.style.backgroundColor = background != null ? background : "";
        classes(this.el).clear().add(...this.css_classes());
        let header_background = this.model.header_background;
        if (!this.model.collapsed && this.model.active_header_background != null)
            header_background = this.model.active_header_background;
        const header = this.child_views[0];
        let header_el;
        if (this.model.collapsible) {
            this.button_el = DOM.createElement("button", { type: "button", class: header_css_classes });
            this.button_el.style.backgroundColor = header_background != null ? header_background : "";
            this.button_el.appendChild(header.el);
            const icon = DOM.createElement("p", { class: button_css_classes });
            icon.innerHTML = this.model.collapsed ? "+" : "\u2212";
            this.button_el.appendChild(icon);
            this.button_el.onclick = () => this._toggle_button();
            header_el = this.button_el;
        }
        else {
            header_el = DOM.createElement(header_tag, { class: header_css_classes });
            header_el.style.backgroundColor = header_background != null ? header_background : "";
            header_el.appendChild(header.el);
        }
        header_el.style.color = header_color != null ? header_color : "";
        this.el.appendChild(header_el);
        header.render();
        for (const child_view of this.child_views.slice(1)) {
            if (!this.model.collapsed)
                this.el.appendChild(child_view.el);
            child_view.render();
        }
    }
    _toggle_button() {
        this.model.collapsed = !this.model.collapsed;
    }
    _collapse() {
        this.invalidate_render();
        this.resize_layout();
    }
    _createElement() {
        return DOM.createElement(this.model.tag, { class: this.css_classes() });
    }
}
CardView.__name__ = "CardView";
export class Card extends Column {
    constructor(attrs) {
        super(attrs);
    }
    static init_Card() {
        this.prototype.default_view = CardView;
        this.define({
            active_header_background: [p.String, null],
            button_css_classes: [p.Array, []],
            collapsed: [p.Boolean, true],
            collapsible: [p.Boolean, true],
            header_background: [p.String, null],
            header_color: [p.String, null],
            header_css_classes: [p.Array, []],
            header_tag: [p.String, "div"],
            tag: [p.String, "div"],
        });
    }
}
Card.__name__ = "Card";
Card.__module__ = "panel.models.layout";
Card.init_Card();
//# sourceMappingURL=card.js.map