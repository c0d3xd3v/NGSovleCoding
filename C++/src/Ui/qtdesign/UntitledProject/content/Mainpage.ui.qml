import QtQuick 6.6
import QtQuick.Controls 6.6
import QtQuick.Layouts 6.6

Item {
    visible: true

    BackgroundToolbar {
        id: toolBar
        height: 48
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: parent.top
        anchors.leftMargin: 0
        anchors.rightMargin: 0
        anchors.topMargin: 0
        Layout.fillWidth: true

        RowLayout {
            id: rowLayout
            anchors.fill: parent

            ToolButton {
                id: toolButton1
                width: 42
                height: 42
                text: qsTr("Tool Button")
                topPadding: 1
                bottomPadding: 0
                rightPadding: 1
                leftPadding: 1
                padding: 0
                spacing: 0
                icon.height: 36
                icon.width: 36
                icon.color: "#00eaeaea"
                icon.cache: false
                display: AbstractButton.IconOnly
                icon.source: "/home/kai/Development/github/NGSovleCoding/C++/src/Ui/qtdesign/UntitledProject/content/icons/openfile.svg"
            }

            TestBTN {
                display: AbstractButton.IconOnly
                pageSource: "Renderingpage.ui.qml"
                pageIndex: 1
            }

            TestBTN {
                display: AbstractButton.IconOnly
                pageSource: "Exportpage.ui.qml"
                pageIndex: 1
            }

            Item {
                id: item1
                width: 200
                height: 10
                Layout.fillHeight: true
                Layout.fillWidth: true
            }
        }
    }

    ScrollView {
        id: scrollView
        anchors.left: parent.left
        anchors.right: parent.right
        anchors.top: toolBar.bottom
        anchors.bottom: parent.bottom
        anchors.leftMargin: 0
        anchors.rightMargin: 0
        anchors.topMargin: 0
        anchors.bottomMargin: 0
        contentWidth: scrollView.width

        ColumnLayout {
            id: columnLayout1
            x: 0
            anchors.left: parent.left
            anchors.right: parent.right
            anchors.top: parent.top
            anchors.bottom: parent.bottom
            anchors.leftMargin: 20
            anchors.rightMargin: 20
            anchors.topMargin: 20
            anchors.bottomMargin: 20
            spacing: 3

            MeshingPane {
                id: columnLayout
                Layout.fillHeight: true
                Layout.fillWidth: true
            }
        }
    }
}
