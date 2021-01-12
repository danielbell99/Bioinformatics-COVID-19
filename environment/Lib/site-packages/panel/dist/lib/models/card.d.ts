import { Column, ColumnView } from "@bokehjs/models/layouts/column";
import { Layoutable } from "@bokehjs/core/layout/layoutable";
import { Column as ColumnLayout } from "@bokehjs/core/layout/grid";
import { Size } from "@bokehjs/core/layout/types";
import * as p from "@bokehjs/core/properties";
export declare class CollapseableColumnLayout extends ColumnLayout {
    collapsed: boolean;
    constructor(items: Layoutable[], collapsed?: boolean);
    protected _measure_totals(row_heights: number[], col_widths: number[]): Size;
}
export declare class CardView extends ColumnView {
    model: Card;
    button_el: HTMLButtonElement;
    connect_signals(): void;
    _update_layout(): void;
    render(): void;
    _toggle_button(): void;
    _collapse(): void;
    protected _createElement(): HTMLElement;
}
export declare namespace Card {
    type Attrs = p.AttrsOf<Props>;
    type Props = Column.Props & {
        active_header_background: p.Property<string | null>;
        button_css_classes: p.Property<string[]>;
        collapsed: p.Property<boolean>;
        collapsible: p.Property<boolean>;
        header_background: p.Property<string | null>;
        header_color: p.Property<string | null>;
        header_css_classes: p.Property<string[]>;
        header_tag: p.Property<string>;
        tag: p.Property<string>;
    };
}
export interface Card extends Card.Attrs {
}
export declare class Card extends Column {
    properties: Card.Props;
    constructor(attrs?: Partial<Card.Attrs>);
    static __module__: string;
    static init_Card(): void;
}
